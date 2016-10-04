
"""
Classes to deal with homogenisation models written in Stan.
"""

import cPickle as pickle
import logging
import numpy as np
import os
import pystan as stan
from astropy.table import Table
from itertools import combinations

from . import plot

logger = logging.getLogger("ges")


def _guess_parameter_name(table, name):
    if name in table.dtype.names:
        p_name = name
    else:
        p_name = name.upper() if name == name.lower() else name.lower()

    guesses = ("e_{}", "E_{}", "{}_ERR", "{}_err", "{}_ERROR", "{}_error")
    for guess in guesses:
        if guess.format(p_name) in table.dtype.names:
            p_e_name = guess.format(p_name)
            break

    else:
        raise ValueError(
            "cannot guess column name for the error in {}".format(name))

    return (p_name, p_e_name)


def _guess_parameter_names(table, names):
    return zip(*map(lambda name: _guess_parameter_name(table, name), names))



def _homogenise_survey_measurements(database, wg, parameter, cname, N=100,
    stan_model=None, stan_chains=None, stan_data=None, update_database=True):
    """
    Produce an unbiased estimate of an astrophyiscal parameter for a given
    survey object.

    :param cname:
        The CNAME (unique star identifier) of an object.

    :param wg:
        The working group to consider measurements from.

    :param parameter:
        The name of the parameter to estimate.

    :param N: [optional]
        The number of samples to draw. Setting to zero or outside the valid
        range defaults to the maximum number of draws availiable.

    :param stan_model: [optional]
        The fitted Stan model. 

        Either the fitted stan model must be provided, or the Stan chains and 
        the data dictionary (`stan_chains` and `stan_data`) must be provided.

    :param stan_chains: [optional]
        The extracted samples (chains) from the sampled Stan model.

        Either the `stan_model` must be provided, or the Stan chains and
        the data dictionary (`stan_chains` and `stan_data`) must be provided.

    :param stan_data: [optional]
        The data dictionary provided to Stan.

        Either the `stan_model` must be provided, or the Stan chains and
        the data dictionary (`stan_chains` and `stan_data`) must be provided.    
    """

    # Get the data for this object.
    estimates = database.retrieve_table(
        """ SELECT  DISTINCT ON (filename, node_id)
                    results.id, cname, node_id, snr, trim(filename) as filename, 
                    {parameter}
            FROM    results, nodes
            WHERE   nodes.wg = {wg}
              AND   nodes.id = results.node_id
              AND   cname = '{cname}'
              AND   {parameter} <> 'NaN'
              AND   passed_quality_control = true;
        """.format(wg=wg, cname=cname, parameter=parameter))

    if estimates is None:
        return (np.array([np.nan]), np.array([np.nan]))

    # Extract N samples for all the parameters.

    # For each sample, calculate:
    #   1. The total variance (systematic**2 + (alpha/SNR)**2)
    #   2. The weighted mean from all observations by that nodes.
    #    --> check that this follows 1/sqrt(N)
    #   3. Construct a covariance matrix using the weighted means, uncertainties
    #       and the correlation coefficients
    #   4. Draw from a Gaussian using the weighted means and your new Cov matrix
    #   5. Record the draw.

    if stan_model is not None:
        samples = stan_model.extract(
            pars=set(stan_model.model_pars).difference(("Sigma", )))
        unique_node_ids = stan_model.data["node_ids"]

    elif stan_chains is not None and stan_data is not None:
        samples, unique_node_ids = (stan_chains, stan_data["node_ids"])

    else:
        raise ValueError(
            "stan_model must be provided, or stan_chains and stan_data")

    K = len(samples["var_intrinsic"])
    M = len(set(estimates["node_id"]))

    if 1 > N or N > K:
        N = K
        indices = range(K)

    else:
        # Select N random indices from K
        indices = np.random.choice(K, N, replace=False)

    estimates = estimates.group_by("node_id")
        
    mu_samples = np.zeros(N)
    var_samples = np.zeros(N)

    rho_terms = list(combinations(unique_node_ids, 2))

    for ii, i in enumerate(indices):

        # Average from each node first.
        mu_node = np.zeros(M)
        var_node_sys = np.zeros(M)
        var_node_rand = np.zeros(M)

        node_ids = np.zeros(M)

        for j, s in enumerate(estimates.groups.indices[:-1]):
            e = estimates.groups.indices[j + 1]
            L = e - s

            k = np.where(estimates["node_id"][s] == unique_node_ids)[0][0]
            node_ids[j] = estimates["node_id"][s]

            # Calculate weighted node mean from random uncertainties.
            mu = np.array(
                estimates[parameter][s:e] + samples["c0_estimators"][i, k])

            spectrum_snr = np.clip(estimates["snr"][s:e], 0.1, np.inf)
            C = np.eye(L) * (samples["alpha_sq"][i, k]/spectrum_snr**2)

            W = np.ones((L, 1))
            Cinv = np.linalg.inv(C)
            var_node_rand[j] = 1.0/np.dot(np.dot(W.T, Cinv), W)
            mu_node[j] = var_node_rand[j] * np.dot(np.dot(W.T, Cinv), mu)

            var_node_sys[j] = samples["var_sys_estimator"][i, k]

        node_ids = node_ids.astype(int)

        # Construct the covariance matrix for node-to-node measurements.
        # (This includes the systematic component.)
        C = np.eye(M) * (var_node_rand + var_node_sys)

        for j in range(M):
            for k in range(j + 1, M):
                # Need the right index to match to rho_estimators.
                l = rho_terms.index(tuple(node_ids[[j, k]]))
                term = samples["rho_estimators"][i, l] * (C[j,j] * C[k,k])**0.5
                C[j, k] = C[k, j] = term

        # Gauss-Markov theorem.
        W = np.ones((M, 1))
        Cinv = np.linalg.inv(C)
        var = 1.0/np.dot(np.dot(W.T, Cinv), W)
        mu = var * np.dot(np.dot(W.T, Cinv), mu_node)

        mu_samples[ii] = mu
        var_samples[ii] = var

    central = np.median(mu_samples)
    error = np.sqrt(np.median(var_samples))

    if update_database:

        # Check if there is an entry for this (wg, cname) in the wg_recommended 
        # table
        #q = np.percentile(posterior, [16, 50, 84])
        #central, pos_error, neg_error = (q[1], q[2] - q[1], q[0] - q[1])

        data = {
            "wg": wg, 
            "cname": cname, 
            "snr": np.nanmedian(estimates["snr"]),
            parameter: central, 
            #"e_pos_{}".format(parameter): pos_error, 
            #"e_neg_{}".format(parameter): neg_error, 
            #"e_{}".format(parameter): np.abs([pos_error, neg_error]).max(),
            "e_{}".format(parameter): error,
            "nn_nodes_{}".format(parameter): len(set(estimates["node_id"])),
            "nn_spectra_{}".format(parameter): len(set(estimates["filename"])),
            "provenance_ids_for_{}".format(parameter): list(estimates["id"].data.astype(int))
        }

        record = database.retrieve(
            """ SELECT id
                  FROM recommended_results
                 WHERE wg = %s
                   AND cname = %s
            """, (wg, cname, ))

        if record:
            database.update(
                """ UPDATE recommended_results
                       SET {}
                     WHERE id = '{}'""".format(
                        ", ".join([" = ".join([k, "%s"]) for k in data.keys()]),
                        record[0][0]),
                data.values())
            logger.info(
                "Updated record {} in recommended_results".format(record[0][0]))

        else:
            new_record = database.retrieve(
                """ INSERT INTO recommended_results ({})
                    VALUES ({}) RETURNING id""".format(
                        ", ".join(data.keys()), ", ".join(["%s"] * len(data))),
                data.values())
            logger.info(
                "Created new record {} in recommended_results ({} / {})"\
                .format(new_record[0][0], wg, cname))

    return (mu_samples, var_samples)







class BaseEnsembleModel(object):

    def __init__(self, database, wg, parameter, calibrators, recompile=False, 
        overwrite=True, **kwargs):
        self._database = database
        self._wg = wg
        self._parameter = parameter
        self._calibrators = calibrators

        self._model = self._load_model(recompile, overwrite, **kwargs)

        # For later.
        self._data = None
        self._fit = None
        self._chains = None

        return None


    def _load_model(self, recompile, overwrite, path=None):
        """
        Load the model.

        :param recompile:
            Re-compile the model if it has already been compiled.

        :param overwrite:
            Overwrite the compiled model path if it already exists.

        :param path: [optional]
            The path of the file containing the model code. Defaults to
            `self._MODEL_PATH`
        """

        # Is the model already compiled?
        path = path or self._MODEL_PATH
        compiled_path = "{}.compiled".format(path)

        while os.path.exists(compiled_path) and not recompile:

            with open(compiled_path, "rb") as fp:
                model = pickle.load(fp)

            # Check that the model code is the same as what we expected.
            with open(path, "r") as fp:
                model_code = fp.read()

            if model_code != model.model_code:
                logger.warn("Pre-compiled model differs to the code in {}; "\
                    "re-compiling".format(path))
                recompile = True
                continue

            else:
                logger.info(
                    "Using pre-compiled model from {}".format(compiled_path))
                break

        else:
            logger.info("Compiling model from {}".format(path))

            model = stan.StanModel(file=path)

            # Save the compiled model.
            if not os.path.exists(compiled_path) or overwrite:
                with open(compiled_path, "wb") as fp:
                    pickle.dump(model, fp, -1)

        return model


    def _validate_stan_inputs(self, **kwargs):
        """
        Check the format of the initial values for the model. If a dictionary
        is specified and multiple chains are given, then the initial value will
        be re-cast as a list of dictionaries (one per chain).
        """

        # Copy the dictionary of keywords.
        kwds = {}
        kwds.update(kwargs)

        # Allow for a keyword that will disable any verification checks.
        if not kwds.pop("validate", True):
            return kwds

        # Check chains and init values.
        if "init" in kwds.keys() and isinstance(kwds["init"], dict) \
        and kwds.get("chains", 1) > 1:

            init, chains = (kwds["init"], kwds.get("chains", 1))
            logger.info(
                "Re-specifying initial values to be list of dictionaries, "\
                "allowing one dictionary per chain ({}). "\
                "Specify validate=False to disable this behaviour"\
                .format(chains))
            
            kwds["init"] = [init] * chains

        if kwargs.get("data", None) is None:
            try:
                self._data 
            except AttributeError:
                self._data, self._metadata = self._prepare_data()

            kwds["data"] = self._data

        return kwds


    def optimize(self, data=None, recompile=False, overwrite=False, **kwargs):
        """
        Optimize the model given the data. Keyword arguments are passed directly
        to the `StanModel.optimizing` method.

        :param data:
            A dictionary containing the required key/value pairs for the STAN
            model.

        :param recompile: [optional]
            Re-compile the model if it has already been compiled.

        :param overwrite: [optional]
            Overwrite the compiled model path if it already exists.
        """

        kwds = self._validate_stan_inputs(data=data, **kwargs)
        return self._model.optimizing(**kwds)


    def sample(self, data=None, chains=4, iter=2000, warmup=None, recompile=False, 
        overwrite=False, **kwargs):
        """
        Draw samples from the model. Keyword arguments are passed directly to
        `StanModel.sampling`.

        :param data:
            A dictionary containing the required key/value pairs for the Stan
            model.

        :param chains: [optional]
            Positive integer specifying the number of chains.

        :param iter:
            Positive integer specifying how many iterations for each chain
            including warmup.

        :param warmup: [optional]
            Positive integer specifying the number of warmup (aka burn-in)
            iterations. As warm-up also specifies the number of iterations used
            for step-size adaption, warmup samples should not be used for
            inference. Defaults to iter // 2.

        :param recompile: [optional]
            Re-compile the model if it has already been compiled.

        :param overwrite: [optional]
            Overwrite the compiled model path if it already exists.
        """

        kwds = self._validate_stan_inputs(
            data=data, chains=chains, iter=iter, warmup=warmup, **kwargs)

        self._fit = self._model.sampling(**kwds)
        return self._fit


    def homogenise_star(self, cname, **kwargs):
        """
        Produce an unbiased estimate of an astrophyiscal parameter for a given
        survey object.

        :param cname:
            The CNAME (unique star identifier) of an object.

        :returns:
            An array of draws from the marginalized posterior distribution.
        """

        return _homogenise_survey_measurements(
            self._database, self._wg, self._parameter, cname,
            stan_chains=self._chains, stan_data=self._data,
            **kwargs)


    def homogenise_all_stars(self, **kwargs):
        """
        Homogenise the stellar astrophysical parameter for all stars in the
        database that are analysed by the current working group. 

        Note: this will homogenise a single parameter for each survey object.
        """

        # Get all unique cnames.
        records = self._database.retrieve_table(
            """ WITH s AS (
                    SELECT id FROM nodes WHERE wg = %s)
                SELECT DISTINCT ON (r.cname) r.cname
                FROM   s, results AS r
            WHERE r.node_id = s.id
            ORDER BY cname DESC""", (self._wg, ))

        # Get samples and data dictionary -- it will be faster.
        if self._chains is None:
            ignore_model_pars = kwargs.get("__ignore_model_pars", ("Sigma", ))
            model_pars = set(self._fit.model_pars).difference(ignore_model_pars)
            self._chains = self._fit.extract(pars=model_pars)

        assert self._data is not None

        N = len(records)
        for i, cname in enumerate(records["cname"]):

            mu, var = _homogenise_survey_measurements(
                self._database, self._wg, self._parameter, cname,
                stan_chains=self._chains, stan_data=self._data,
                **kwargs)
            pos_error = neg_error = np.median(var**0.5)

            #p = np.percentile(posterior, [16, 50, 84])
            #value, pos_error, neg_error = (p[1], p[2] - p[1], p[0] - p[1])

            logger.info("Homogenised {parameter} for {cname} (WG{wg} {i}/{N}): "
                "{mu:.2f} ({pos_error:.2f}, {neg_error:.2f})".format(
                    parameter=self._parameter, cname=cname, wg=self._wg, i=i+1,
                    N=N, mu=np.median(mu), pos_error=pos_error, neg_error=neg_error))

        if kwargs.get("update_database", True):
            self._database.connection.commit()

        return None


    def write(self, filename, overwrite=False, **kwargs):
        """
        Write the model to disk, including any MCMC chains and data dictionaries
        from Stan.

        :param filename:
            The local path of where to write the model.

        :param overwrite: [optional]
            Overwrite the `filename` if it already exists.
        """

        if os.path.exists(filename) and not overwrite:
            raise IOError(
                "filename {} exists and not overwriting".format(filename))

        state = {
            "model_path": self._MODEL_PATH, 
            "wg": self._wg,
            "parameter": self._parameter,
            "calibrators": self._calibrators,
            "data": self._data,
            "metadata": self._metadata,
            "chains": None
        }

        if self._chains is not None:
            state["chains"] = self._chains

        elif self._fit is not None:
            # Extract chains.
            ignore_model_pars = kwargs.get("__ignore_model_pars", ("Sigma", ))
            model_pars = set(self._fit.model_pars).difference(ignore_model_pars)
            self._chains = self._fit.extract(pars=model_pars)

            state["chains"] = self._chains

        with open(filename, "wb") as fp:
            pickle.dump(state, fp, -1)

        return None


    @classmethod
    def read(cls, filename, database, **kwargs):
        """
        Read a saved model from disk.

        :param filename:
            The local path of where the model is saved.

        :param database:
            The database connection that this model will use.
        """

        with open(filename, "rb") as fp:
            state = pickle.load(fp)

        klass = \
            cls(database, state["wg"], state["parameter"], state["calibrators"],
                **kwargs)

        # Update the klass.
        klass._data = state.get("data", None)
        klass._metadata = state.get("metadata", None)
        klass._chains = state.get("chains", None)

        return klass


    def node_uncertainty_with_snr(self, **kwargs):
        return plot.node_uncertainty_with_snr(self, **kwargs)




class MeanEnsembleModel(BaseEnsembleModel):

    _MODEL_PATH = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "mean-ensemble-model-with-correlations.stan")

    def _prepare_data(self, parameter=None, default_sigma_calibrator=1e3,
        fill_function=np.mean, fill_variance=1e50, require_no_gaps=False,
        minimum_snr=1, maximum_snr=500, 
        sql_constraint=None):
        """
        Prepare the data for the model so that it can be supplied to Stan.

        """        
        parameter = str(parameter or self._parameter).lower()
        if parameter not in ("teff", "logg", "feh", "xi", "mh"):
            raise AreYouSureYouKnowWhatYoureDoing

        # Get the data from the database for this WG.
        # TODO: Need a better way to identify calibrators. Right now we do it
        #       just on GES_TYPE, but in future we may want to do this directly
        #       on to CNAMEs in the calibrator list.
        data = self._database.retrieve_table(
            """ WITH    n AS (
                            SELECT id FROM nodes WHERE wg = {wg}), 
                        s AS (
                            SELECT cname, ges_type, ges_fld
                            FROM spectra
                            WHERE ges_type LIKE 'GE_SD_B%') 
                SELECT DISTINCT ON (r.filename, r.node_id) 
                    s.cname, s.ges_type, s.ges_fld,
                    r.filename, r.node_id, r.snr,
                    r.{parameter}, r.e_{parameter}
                FROM    s, n, results as r 
                WHERE r.cname = s.cname 
                  AND r.node_id = n.id 
                  AND r.passed_quality_control = true {sql_constraint}
                """.format(
                    wg=self._wg, parameter=parameter,
                    sql_constraint=""   if sql_constraint is None \
                                        else " AND {}".format(sql_constraint)))
        assert data is not None, "No calibrator data from WG {}".format(wg)

        # Calibrator parameter name
        calibrator_name, calibrator_e_name = _guess_parameter_name(
            self._calibrators, parameter)

        # Match calibrators to data.
        data_ges_flds = map(str.strip, data["ges_fld"])
        common_ges_flds \
            = set(data_ges_flds).intersection(self._calibrators["GES_FLD"])
        keep = np.array([each in common_ges_flds for each in data_ges_flds])
        data = data[keep]

        unique_cnames = np.sort(np.unique(data["cname"]))
        unique_node_ids = np.sort(np.unique(data["node_id"]))
        calibrator_indices = np.nan * np.ones(len(unique_cnames))

        for i, cname in enumerate(unique_cnames):

            ges_fld = data["ges_fld"][data["cname"] == cname][0]
            match = (ges_fld == self._calibrators["GES_FLD"][i])
            calibrator_indices[i] = np.where(match)[0][0]

        calibrators = self._calibrators[calibrator_indices]
        
        data = data.group_by(["cname", "node_id"])

        # Get the maximum number of visits per calibrator.
        max_visits = max(map(len, data.groups))

        shape = (N_calibrators, N_nodes, max_visits)
        estimates = np.zeros(shape)
        ivar_spectrum = np.zeros(shape)

        N_visits = np.zeros((N_calibrators, N_nodes))

        for i, si in enumerate(data.groups.indices[:-1]):
            ei = data.groups.indices[i + 1]

            j = np.where(data["cname"][si] == unique_cnames)[0][0]
            k = np.where(data["node_id"][si] == unique_node_ids)[0][0]

            N = ei - si
            N_visits[j, k] = N

            snr = np.clip(data["snr"][si:ei], 1, np.inf)
            ivar_spectrum[j, k, :N] = 1.0/snr**2

            estimates[j, k, :N] = data[parameter][si:ei]


        N_calibrators, N_nodes = len(unique_cnames), len(unique_node_ids)
        data_dict = {
            "N_calibrators": N_calibrators,
            "N_nodes": N_nodes,
            "N_visits": N_visits,
            "max_visits": max_visits,
            
            "mu_calibrator": np.array(calibrators[calibrator_name]),
            "sigma_calibrator": np.array(calibrators[calibrator_e_name]),


            "ivar_spectrum": ivar_spectrum,
            "estimates": estimates,

            # Include the node_ids so that we will absolutely have them when we
            # try to run the homogenisation.
            "node_ids": unique_node_ids,
            "calibrator_cnames": unique_cnames,
        }
        if not np.all(np.isfinite(data_dict["sigma_calibrator"])):
            logger.warn(
                "Not all finite calibrator entries of {} have finite errors! " \
                "Setting default sigma_calibrator as {}".format(
                    parameter, default_sigma_calibrator))
            finite = np.isfinite(data_dict["sigma_calibrator"])
            data_dict["sigma_calibrator"][~finite] = default_sigma_calibrator


        # Prepare metadata.
        node_names = self._database.retrieve_table("SELECT * FROM nodes")
        metadata = {
            "N_pairwise_nodes": np.math.factorial(N_nodes) \
                / (2*np.math.factorial(max(1, N_nodes - 2))),
            "calibrator_ges_fld": calibrators["ges_fld"],
            "node_ids": unique_node_ids,
            "node_names": \
                [node_names["name"][node_names["id"] == node_id][0].strip() \
                    for node_id in unique_node_ids]
        }
        result = (data_dict, metadata)

        self._data, self._metadata = result
        return result



class SingleParameterEnsembleModel(BaseEnsembleModel):

    _MODEL_PATH = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "single-parameter-ensemble-model-with-correlations.stan")

    def _prepare_data(self, parameter=None, default_sigma_calibrator=1e3,
        fill_function=np.mean, fill_variance=1e50, require_no_gaps=False,
        minimum_snr=1, maximum_snr=500, include_calibrator_function=None, 
        sql_constraint=None):
        """
        Prepare the data for the model so that it can be supplied to Stan.

        :param parameter: [optional]
            The name of the model parameter (e.g., teff) that will be used in
            this single parameter ensemble model. If `None` provided, then this
            defaults to the model parameter provided when the EnsembleModel was
            initiated.

        :param fill_function: [optional]
            A function that takes as input finite values and returns as output a
            value that will be used to fill arrays that would otherwise be non-
            finite entries. As long as the `fill_variance` is suitably large, it
            does not matter what the fill value is.

        :param fill_variance: [optional]
            The variance to add for matrix entries that have been filled with
            the `fill_function` value.

        :param require_no_gaps: [optional]
            Only include observations where there is a measurement from every
            estimator (node). This should *only* be used for testing models.

        :param include_calibrator_function: [optional]
            A function that takes as input a row from a table containing the
            information about the calibrator, and returns a boolean if that
            calibrator should be included in the data for the model. This can
            be used to exclude specific calibrators during cross-validation
            tests.

        :param sql_constraint: [optional]
            Specify an optional SQL constraint to apply.
        """

        parameter = str(parameter or self._parameter).lower()
        if parameter not in ("teff", "logg", "feh", "xi", "mh"):
            raise AreYouSureYouKnowWhatYoureDoing

        # Get the data from the database for this WG.
        # TODO: Need a better way to identify calibrators. Right now we do it
        #       just on GES_TYPE, but in future we may want to do this directly
        #       on to CNAMEs in the calibrator list.
        data = self._database.retrieve_table(
            """ WITH    n AS (
                            SELECT id FROM nodes WHERE wg = {wg}), 
                        s AS (
                            SELECT cname, ges_type, ges_fld
                            FROM spectra
                            WHERE ges_type LIKE 'GE_SD_B%') 
                SELECT DISTINCT ON (r.filename, r.node_id) 
                    s.cname, s.ges_type, s.ges_fld,
                    r.filename, r.node_id, r.snr,
                    r.{parameter}, r.e_{parameter}
                FROM    s, n, results as r 
                WHERE r.cname = s.cname 
                  AND r.node_id = n.id 
                  AND r.passed_quality_control = true {sql_constraint}
                """.format(
                    wg=self._wg, parameter=parameter,
                    sql_constraint=""   if sql_constraint is None \
                                        else " AND {}".format(sql_constraint)))
        assert data is not None, "No calibrator data from WG {}".format(wg)

        # Calibrator parameter name
        calibrator_name, calibrator_e_name = _guess_parameter_name(
            self._calibrators, parameter)

        finite_calibrator = np.isfinite(self._calibrators[calibrator_name])
        if not np.all(finite_calibrator):
            logger.warn("Not all calibrator values of {} are finite! ({}/{})"\
                .format(parameter, sum(finite_calibrator), len(finite_calibrator)))

        calibrators = self._calibrators[finite_calibrator]
        if include_calibrator_function is not None:
            keep = np.array([include_calibrator_function(row) for row in calibrators])
            logger.info(
                "Excluded {} calibrators based on include_calibrator_function"\
                .format(len(keep) - sum(keep)))
            calibrators = calibrators[keep]

        N_calibrators = len(calibrators)
        N_estimators = len(set(data["node_id"]))
        data_dict = {
            "N_calibrators": N_calibrators,
            "N_estimators": N_estimators,
            "mu_calibrator": np.array(calibrators[calibrator_name]),
            "sigma_calibrator": np.array(calibrators[calibrator_e_name]),
        }

        # Check the data so far.
        if 1 > N_calibrators:
            raise ValueError("no calibrators")

        if 1 > N_estimators:
            raise ValueError("no estimators")

        if not np.all(np.isfinite(data_dict["sigma_calibrator"])):
            logger.warn(
                "Not all finite calibrator entries of {} have finite errors! " \
                "Setting default sigma_calibrator as {}".format(
                    parameter, default_sigma_calibrator))
            finite = np.isfinite(data_dict["sigma_calibrator"])
            data_dict["sigma_calibrator"][~finite] = default_sigma_calibrator

        # OK now update the data dictionary with the spectroscopic measurements.
        # Need to group by node id and CNAME.
        unique_calibrators = np.array(map(str.strip, calibrators["GES_FLD"]))
        unique_estimators = np.sort(np.unique(data["node_id"]))

        assert len(set(calibrators["GES_FLD"])) == N_calibrators
        assert len(unique_estimators) == N_estimators

        # Group data by FILENAME.
        skipped = {}
        data = data.group_by(["filename"])
        N_calibrator_visits = len(data.groups)

        spectrum_snr = np.nan * np.ones(N_calibrator_visits)
        calibrator_index = -1 * np.ones(N_calibrator_visits, dtype=int)
        estimates = np.nan * np.ones((N_calibrator_visits, N_estimators))

        for i, si in enumerate(data.groups.indices[:-1]):
            ei = data.groups.indices[i + 1]

            ges_fld = data["ges_fld"][si].strip()

            try:
                j = np.where(ges_fld == unique_calibrators)[0][0]

            except IndexError:
                if ges_fld not in skipped:
                    logger.warn("Skipping calibrator {}".format(ges_fld))
                    skipped.setdefault(ges_fld, [])
                skipped[ges_fld].append(i)
                continue

            calibrator_index[i] = j
            spectrum_snr[i] = np.nanmedian(data["snr"][si:ei])
            if np.std(data["snr"][si:ei]) > 1:
                logger.warn(
                    "Standard deviation in SNR for a single {} observation is {:.1f}: {}"\
                    .format(ges_fld, np.std(data["snr"][si:ei]), data["snr"][si:ei]))

            for k in range(si, ei):
                l = np.where(data["node_id"][k] == unique_estimators)[0][0]
                estimates[i, l] = data[parameter][k]
        
        # Fill in the NaN estimates with the mean of that entry, and create an
        # additive variance array.
        if require_no_gaps:
            var_additive = np.zeros_like(estimates)

        else:
            non_finite = ~np.isfinite(estimates)
            var_additive = np.zeros_like(estimates)
            var_additive[non_finite] = fill_variance
            for i in range(estimates.shape[0]):
                fill_value = fill_function(
                    estimates[i, (~non_finite)[i]])
                assert np.isfinite(fill_value) or not any(np.isfinite(estimates[i]))
                estimates[i, non_finite[i]] = fill_value

        # Only keep entries where we have *any* finite estimates.
        # This should also remove any skipped calibrators.
        keep = np.all(np.isfinite(estimates), axis=1) \
            * (spectrum_snr >= minimum_snr) \
            * (maximum_snr >= spectrum_snr)
        
        # Ensure spectrum_snr values are valid.
        #spectrum_snr = np.clip(spectrum_snr, minimum_snr, maximum_snr)
        data_dict.update({
            "ivar_spectrum": (1.0/spectrum_snr[keep])**2,
            "N_calibrator_visits": sum(keep),
            "var_additive": var_additive[keep],
            "estimates": estimates[keep],
            # Make the calibrator index 1-indexed for STAN
            "calibrator_index": 1 + calibrator_index[keep].astype(int),

            # Include the node_ids so that we will absolutely have them when we
            # try to run the homogenisation.
            "node_ids": unique_estimators,
        })

        Sigma_c0 = np.array(map(np.diag, var_additive[keep]))
        Sigma_c1 = np.ones_like(Sigma_c0)
        for i in range(sum(keep)):
            for j in range(N_estimators):
                if Sigma_c0[i, j, j] > 0:
                    Sigma_c1[i, :, j] = 0
                    Sigma_c1[i, j, :] = 0

        data_dict.update({
            "Sigma_c0": Sigma_c0,
            "Sigma_c1": Sigma_c1
        })

        # Create additional metadata
        node_names = self._database.retrieve_table("SELECT * FROM nodes")
        metadata = {
            "N_pairwise_estimators": np.math.factorial(N_estimators) \
                / (2*np.math.factorial(max(1, N_estimators - 2))),
            "calibrator_ges_fld": unique_calibrators,
            "node_ids": unique_estimators,
            "node_names": \
                [node_names["name"][node_names["id"] == node_id][0].strip() \
                    for node_id in unique_estimators]
        }
        result = (data_dict, metadata)
        self._data, self._metadata = result
        return result

