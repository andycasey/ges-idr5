
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

    K = samples["truths"].shape[0]
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

        # Weighted mean from each node first (over many observations)
        mu_node = np.zeros(M)
        var_node_sys = np.zeros(M)
        var_node_rand = np.zeros(M)

        node_ids = np.zeros(M)

        for j, s in enumerate(estimates.groups.indices[:-1]):
            e = estimates.groups.indices[j + 1]
            L = e - s

            k = np.where(estimates["node_id"][s] == unique_node_ids)[0][0]
            node_ids[j] = estimates["node_id"][s]

            # Get the 'unbiased' values.
            mu = np.array(
                estimates[parameter][s:e] + samples["c0_estimators"][i, k])

            spectrum_snr = np.clip(estimates["snr"][s:e], 1, np.inf)
            C = np.eye(L) * (samples["alpha"][i, k]**2/spectrum_snr)

            W = np.ones((L, 1))
            Cinv = np.linalg.inv(C)
            
            raise a
            # this is wrong here because we want to just sum the inv vars
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

        if np.linalg.cond(C) > np.finfo(matrix.dtype).eps:
            # Matrix is ill-conditioned.
            logger.warn("Covariance matrix is ill-conditioned.")

            # Use SVD to invert matrix.
            U, s, V = np.linalg.svd(C)
            Cinv = np.dot(np.dot(V.T, np.linalg.inv(np.diag(s))), U.T)

        else:
            Cinv = np.linalg.inv(C)

        var = 1.0/np.dot(np.dot(W.T, Cinv), W)
        if var < 0:
            logger.warn("Negative variance returned!")
            
        mu = np.abs(var) * np.dot(np.dot(W.T, Cinv), mu_node)

        mu_samples[ii] = mu
        var_samples[ii] = var

    central = np.nanmedian(mu_samples)
    error = np.sqrt(np.nanmedian(np.abs(var_samples)))

    if np.isfinite(central):
        assert np.isfinite(error)

    raise a

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

        model_code = kwargs.get("model_code", None)
        if model_code is None:
            with open(self._MODEL_PATH, "r") as fp:
                model_code = fp.read()

        self._model_code = model_code
        
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
            pos_error = neg_error = np.nanmedian(var**0.5)

            #p = np.percentile(posterior, [16, 50, 84])
            #value, pos_error, neg_error = (p[1], p[2] - p[1], p[0] - p[1])

            logger.info("Homogenised {parameter} for {cname} (WG{wg} {i}/{N}): "
                "{mu:.2f} ({pos_error:.2f}, {neg_error:.2f})".format(
                    parameter=self._parameter, cname=cname, wg=self._wg, i=i+1,
                    N=N, mu=np.nanmedian(mu), pos_error=pos_error, neg_error=neg_error))

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
            "model_code": self._model.model_code, 
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
                model_code=state.get("model_code", None), **kwargs)

        # Update the klass.
        klass._data = state.get("data", None)
        klass._metadata = state.get("metadata", None)
        klass._chains = state.get("chains", None)

        return klass


    def plot_node_uncertainty_with_snr(self, **kwargs):
        return plot.node_uncertainty_with_snr(self, **kwargs)


    def plot_node_correlations(self, **kwargs):
        return plot.node_correlations(self, **kwargs)



class EnsembleModel(BaseEnsembleModel):

    _MODEL_PATH = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "ensemble-model.stan")

    def _prepare_data(self, parameter=None, default_sigma_calibrator=1e3,
        minimum_node_estimates=1, sql_constraint=None):
        """
        Prepare the data for the model so that it can be supplied to Stan.

        :param parameter: [optional]
            The name of the model parameter (e.g., teff) that will be used in
            this single parameter ensemble model. If `None` provided, then this
            defaults to the model parameter provided when the EnsembleModel was
            initiated.

        :param minimum_node_estimates: [optional]
            The minimum number of node measurements for a single visit of a 
            calibrator spectrum. If this is set to a negative value, then only
            calibrator results will be used where *all* nodes have provided a
            measurement for that spectrum.

        :param sql_constraint: [optional]
            Specify an optional SQL constraint to include when retrieving the
            node results.
        """

        parameter = str(parameter or self._parameter).lower()
        valid_parameters = ["teff", "logg", "feh"]
        if parameter not in valid_parameters:
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

        # Calibrator parameter names
        calibrator_name, calibrator_e_name = _guess_parameter_name(
            self._calibrators, parameter)

        finite_calibrator = np.isfinite(self._calibrators[calibrator_name])
        if not np.all(finite_calibrator):
            logger.warn("Not all calibrator values of {} are finite! ({}/{})"\
                .format(parameter, sum(finite_calibrator), len(finite_calibrator)))
       
        # OK now update the data dictionary with the spectroscopic measurements.
        # Need to group by node id and CNAME.
        calibrators = self._calibrators.copy()

        # Common calibrators to serve as an indexing reference.
        common_calibrators = set(map(str.strip, calibrators["GES_FLD"]))\
                                .intersection(map(str.strip, data["ges_fld"]))
        common_calibrators = np.sort(list(common_calibrators))

        # Remove calibrators not in common.
        keep = np.ones(len(self._calibrators), dtype=bool)
        for i, ges_fld in enumerate(self._calibrators["GES_FLD"]):
            if ges_fld.strip() not in common_calibrators:
                keep[i] = False

        calibrators = self._calibrators[keep]
        calibrators.sort("GES_FLD")
        assert calibrators["GES_FLD"][0].strip() == common_calibrators[0]
        C = len(calibrators)

        unique_estimators = np.sort(np.unique(data["node_id"]))
        N = unique_estimators.size

        # Get the maximum number of visits for any calibrator
        data = data.group_by(["cname"])
        V = np.max([len(set(group["filename"])) for group in data.groups])

        skipped = {}
        visits = np.zeros(C, dtype=int)
        estimates = np.nan * np.ones((C, N, V))
        spectrum_ivar = np.zeros((C, V))

        for i, si in enumerate(data.groups.indices[:-1]):
            ges_fld = data["ges_fld"][si].strip()
            if ges_fld not in common_calibrators: continue

            c = np.where(ges_fld == common_calibrators)[0][0]

            ei = data.groups.indices[i + 1]
            filename_visits = np.unique(data["filename"][si:ei])

            for k in range(si, ei):

                n = np.where(data["node_id"][k] == unique_estimators)[0][0]
                v = np.where(data["filename"][k] == filename_visits)[0][0]

                estimates[c, n, v] = data[parameter][k]

                snr = np.clip(data["snr"][k], 1, np.inf)
                spectrum_ivar[c, v] = snr**2

            visits[c] = np.sum(np.any(np.isfinite(estimates[c]), axis=0))

        
        # Remove any nodes with zero measurements.
        keep = np.array([np.any(np.isfinite(estimates[:, n, :])) for n in range(N)])
        estimates = estimates[:, keep, :]
        unique_estimators = unique_estimators[keep]
        N = sum(keep)

        if minimum_node_estimates != 0:
            if minimum_node_estimates < 0:
                minimum_node_estimates = N

            # Only keep visits that have at least the number of minimum node measurements
            for c in range(C):
                mask = np.sum(np.isfinite(estimates[c]), axis=0) >= minimum_node_estimates
                
                n_full_rank = mask.sum()
                estimates[c][:, :n_full_rank] = estimates[c][:, mask]
                estimates[c][:, n_full_rank:] = np.nan

                spectrum_ivar[c, :n_full_rank] = spectrum_ivar[c][mask]
                spectrum_ivar[c, n_full_rank:] = 0

                # Update the number of visits for this calibrator
                visits[c] = n_full_rank

            _slice = visits > 0

            calibrators = calibrators[_slice]

            visits = visits[_slice]
            estimates = estimates[_slice, :, :max(visits)]
            spectrum_ivar = spectrum_ivar[_slice, :max(visits)]

            C, _, V = estimates.shape

        # Construct N_missing
        is_missing = np.zeros(estimates.shape, dtype=bool)
        for c, v in enumerate(visits):
            is_missing[c] = (~np.isfinite(estimates[c])) \
                          * (np.tile(np.arange(V), N).reshape(N, V) < v)

        if minimum_node_estimates >= N: assert np.sum(is_missing) == 0

        mu_calibrator = np.array(calibrators[calibrator_name])
        sigma_calibrator = np.array(calibrators[calibrator_e_name])

        if not np.all(np.isfinite(sigma_calibrator)):
            logger.warn("Not all calibrator uncertainties are finite! "
                        "Filling in with default value {}".format(
                            default_sigma_calibrator))
            sigma_calibrator[~np.isfinite(sigma_calibrator)] = default_sigma_calibrator

        data_dict = {
            "N": N, # number of nodes
            "C": C, # number of calibrators
            "V": V, # maximum number of visits to any calibrator
            "visits": visits.astype(int),

            "is_missing": is_missing.astype(int),
            "TM": np.sum(is_missing),

            "estimates": estimates,
            "spectrum_ivar": spectrum_ivar.T,
            "spectrum_isnr": 1.0/np.sqrt(spectrum_ivar.T),

            "mu_calibrator": mu_calibrator,
            "sigma_calibrator": sigma_calibrator,

            "S": 3, # TODO: Make this flexible? Make the calibrator values accessible from kwargs?
            "all_mu_calibrator": np.vstack(
                [calibrators[p] for p in ("TEFF", "LOGG", "FEH")]).T
        }

        alpha_bounds = dict(teff=(100, 1000), logg=(0.1, 1.0), feh=(0.1, 1.0))
        data_dict.update(
            dict(zip(("lower_alpha_sq", "upper_alpha_sq"), np.array(alpha_bounds[parameter])**2)))

        bounds = dict(teff=(3000, 8000), logg=(0, 5), feh=(-3.5, 0.5))
        data_dict.update(
            dict(zip(("lower_bound", "upper_bound"), bounds[parameter])))

        # Create additional metadata
        node_names = self._database.retrieve_table("SELECT * FROM nodes")
        metadata = {
            "calibrators": calibrators,
            "node_ids": unique_estimators,
            "node_names": \
                [node_names["name"][node_names["id"] == node_id][0].strip() \
                    for node_id in unique_estimators]
        }
        result = (data_dict, metadata)
        self._data, self._metadata = result
        return result