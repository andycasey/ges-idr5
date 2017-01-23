
"""
Classes to perform homogenisation using both a pre-mapping of parameters onto a
`reference node', and a noise model.
"""


import cPickle as pickle
import logging
import numpy as np
import os
import pystan as stan
import scipy.optimize as op
from astropy.table import Table, row

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



class BaseCombinedModel(object):

    def __init__(self, database, wg, parameter, calibrators, recompile=False,
        overwrite=True, **kwargs):
        
        # Store things first.
        self._database = database
        self._wg = wg
        self._parameter = parameter
        self._calibrators = calibrators

        model_code = kwargs.get("model_code", None)
        if model_code is None:
            model_path = kwargs.get("model_path", self._MODEL_PATH)
            with open(model_path, "r") as fp:
                model_code = fp.read()
            self._MODEL_PATH = model_path

        self._model_code = model_code

        # For later.
        self._data = None
        self._fit = None
        self._chains = None

        return None



    def _load_model(self, recompile=False, overwrite=False, path=None, **kwargs):
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

        if getattr(self, "_model", None) is not None:
            return self._model

        # Is the model already compiled?
        path = path or self._MODEL_PATH
        compiled_path = "{}.compiled".format(path)

        while os.path.exists(compiled_path) and not recompile:

            with open(compiled_path, "rb") as fp:
                model = pickle.load(fp)

            # Check that the model code is the same as what we expected.
            assert self._model_code is not None
            
            if self._model_code != model.model_code:
                logger.warn("Pre-compiled model differs to the code in {}; "\
                    "re-compiling".format(path))
                recompile = True
                continue

            else:
                logger.info(
                    "Using pre-compiled model from {}".format(compiled_path))
                break

        else:
            logger.info("Compiling model")

            model = stan.StanModel(model_code=self._model_code)

            # Save the compiled model.
            if not os.path.exists(compiled_path) or overwrite:
                with open(compiled_path, "wb") as fp:
                    pickle.dump(model, fp, -1)

        self._model = model
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

        self._load_model(recompile=recompile, overwrite=overwrite, **kwargs)

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

        self._load_model(recompile=recompile, overwrite=overwrite, **kwargs)

        kwds = self._validate_stan_inputs(
            data=data, chains=chains, iter=iter, warmup=warmup, **kwargs)

        self._fit = self._model.sampling(**kwds)
        return self._fit


    def _extract_chains(self, **kwargs):

        ignore_model_pars = kwargs.get("__ignore_model_pars", ("Sigma", ))
        model_pars = set(self._fit.model_pars).difference(ignore_model_pars)
        self._chains = self._fit.extract(pars=model_pars)

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
            self._extract_chains()
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



class CombinedModel(BaseCombinedModel):

    _MODEL_PATH = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "combined-model.stan")

    
    def _get_mapping_coefficients_from_calibrators(self, vector_terms, 
        sql_constraint=None):
        """
        Map individual node results to the calibrators.

        :param vector_terms:
            A list of two-length tuples containing: (1) the parameter(s) to use
            as a single vector term, and (2) the power for that term.

            For example:

            >>> vector_terms=[
                (("teff", ), 2),
                (("teff", "logg"), 1),
                (("logg", ), 3.5)
            ]

            Would create a vector term with a scaling offset (always present) 
            and terms of `\theta_1 * teff^2`, `\theta_2 * teff * logg`, and
            `\theta_3 * logg^3`.
        
        :param sql_constraint: [optional]
            Specify an optional SQL constraint to include when retrieving the
            node results.
        """

        # What terms should we require be finite?
        require_finite = [self._parameter]
        for terms, power in vector_terms:
            require_finite.extend(list(terms))
        require_finite = list(set(require_finite))

        # Get the node data for all the benchmarks.
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
                        r.teff, r.e_teff,
                        r.logg, r.e_logg,
                        r.feh, r.e_feh
                FROM    s, n, results as r 
                WHERE   r.cname = s.cname 
                  AND   r.node_id = n.id 
                  AND   r.passed_quality_control = true 
                  {finite_constraint}
                  {sql_constraint}
            """.format(wg=self._wg, parameter=self._parameter,
                finite_constraint="".join(
                    ["AND r.{} <> 'NaN' ".format(p) for p in require_finite]),
                sql_constraint="" if sql_constraint is None \
                                  else "AND {}".format(sql_constraint)))
        assert data is not None, "No mapping data from WG{}".format(self._wg)

        # Get the data in a nice format and match to the benchmark value.
        calibrators = self._calibrators.copy()
        calibrator_ges_flds = np.array(map(str.strip, calibrators["GES_FLD"]))

        x, y, before, after = ({}, {}, {}, {})

        data = data.group_by(["ges_fld"])
        G = len(data.groups)

        for i, si in enumerate(data.groups.indices[:-1]):
            ei = data.groups.indices[i + 1]

            print(i, G)

            # Match to the calibrators.
            try:
                k = np.where(calibrator_ges_flds == data["ges_fld"][si].strip())[0][0]

            except IndexError:
                logger.warn("No calibrator values found for {}".format(
                    data["ges_fld"][si]))
                continue # on to the next calibrator group

            else:
                calibrator_value = calibrators[self._parameter.upper()][k]

            for j in range(ei - si):

                # Set up the matrix entries.
                _x = [1]
                for terms, power in vector_terms:
                    value = 1.0
                    for term in terms:
                        value *= data[term][si + j]
                    value **= power
                    _x.append(value)

                node_id = data["node_id"][si + j]
                this_node_value = data[self._parameter][si + j] 

                _y = this_node_value - calibrator_value

                x.setdefault(node_id, [])
                y.setdefault(node_id, [])
                before.setdefault(node_id, [])

                x[node_id].append(_x)
                y[node_id].append(_y)
                before[node_id].append(this_node_value)

        # Prepare the matrices and model.
        for node_id in x.keys():
            x[node_id] = np.atleast_2d(x[node_id]).T
            y[node_id] = np.array(y[node_id])
            before[node_id] = np.array(before[node_id])

        # For each node, find the coefficients to map the results.
        model = lambda x, *theta: np.dot(theta, x)

        mapping_coefficients = {}

        for node_id in set(data["node_id"]):

            xn, yn = (x[node_id], y[node_id])
            theta, cov = op.curve_fit(model, xn, yn, p0=np.zeros(xn.shape[0]))
            
            mapping_coefficients[node_id] = (theta, cov)

            # Calculate post-mapped values.
            after[node_id] = before[node_id] - model(xn, *theta)

            # A sanity check.
            calibrator_values = before[node_id] - yn
            before_stddev = np.std(before[node_id] - calibrator_values)
            after_stddev = np.std(after[node_id] - calibrator_values)

            assert after_stddev <= before_stddev, \
                "You're doing *worse*?! What?!"

            print(node_id, before_stddev, after_stddev)

        return mapping_coefficients


    def _get_mapping_coefficients_from_reference_node(self, vector_terms, 
        reference_node_id, sql_constraint=None):
        """
        Map individual node results to a single reference node.

        :param vector_terms:
            A list of two-length tuples containing: (1) the parameter(s) to use
            as a single vector term, and (2) the power for that term.

            For example:

            >>> vector_terms=[
                (("teff", ), 2),
                (("teff", "logg"), 1),
                (("logg", ), 3.5)
            ]

            Would create a vector term with a scaling offset (always present) 
            and terms of `\theta_1 * teff^2`, `\theta_2 * teff * logg`, and
            `\theta_3 * logg^3`.
        
        :param reference_node_id:
            The id number of the reference node. The reference node must be part
            of the working group for this model.

        :param sql_constraint: [optional]
            Specify an optional SQL constraint to include when retrieving the
            node results.
        """

        assert reference_node_id is not None, "map to benchmarks in this sit."

        # What terms should we require be finite?
        require_finite = [self._parameter]
        for terms, power in vector_terms:
            require_finite.extend(list(terms))
        require_finite = list(set(require_finite))

        data = self._database.retrieve_table(
            """ WITH    n AS (SELECT id FROM nodes WHERE wg = {wg})
                SELECT DISTINCT ON (filename, r.node_id) 
                        TRIM(r.filename) AS filename, r.node_id, r.snr,
                        r.teff, r.e_teff,
                        r.logg, r.e_logg,
                        r.feh, r.e_feh
                FROM    n, results as r 
                WHERE   r.node_id = n.id 
                  AND   r.passed_quality_control = true
                  {sql_constraint}
                  {finite_constraint}
            """.format(wg=self._wg, parameter=self._parameter,
                finite_constraint="".join(
                    ["AND r.{} <> 'NaN' ".format(p) for p in require_finite]),
                sql_constraint="" if sql_constraint is None \
                                  else "AND {}".format(sql_constraint)))        
        assert data is not None, "No mapping data from WG{}".format(self._wg)

        # Now, get the data in a nice format.
        node_ids = set(data["node_id"])
        assert reference_node_id in node_ids, \
            "Reference node with id {} is not in the retrieved results! "\
            "Available nodes are: {}".format(reference_node_id, node_ids)

        # Because we don't know how many results we will have between nodes.
        x, y, before, after = ({}, {}, {}, {})
        node_ids_to_map = node_ids.difference([reference_node_id])
        for node_id in node_ids_to_map:
            x.setdefault(node_id, [])
            y.setdefault(node_id, [])
            before.setdefault(node_id, [])
            after.setdefault(node_id, [])

        # Iterate over the result from each filename.
        data = data.group_by(["filename"])
        G = len(data.groups)

        for i, si in enumerate(data.groups.indices[:-1]):
            ei = data.groups.indices[i + 1]

            print(i, G)

            if 2 > ei - si:
                # Single row -- not useful to us; continue to the next group.
                continue

            try:
                k = np.where(data["node_id"][si:ei] == reference_node_id)[0][0]
            
            except IndexError:
                # If the reference_node_id is not in this group,
                # then we have nothing to compare to.
                continue

            else:
                reference_node_value = data[self._parameter][si + k]


            for j in range(ei - si):
                if j == k: continue
                
                # Set up the matrix entries.
                _x = [1]
                for terms, power in vector_terms:
                    value = 1.0
                    for term in terms:
                        value *= data[term][si + j]
                    value **= power
                    _x.append(value)

                node_id = data["node_id"][si + j]
                this_node_value = data[self._parameter][si + j] 

                _y = this_node_value - reference_node_value

                x[node_id].append(_x)
                y[node_id].append(_y)
                before[node_id].append(this_node_value)


        # Prepare the matrices and model.
        for node_id in x.keys():
            x[node_id] = np.atleast_2d(x[node_id]).T
            y[node_id] = np.array(y[node_id])
            before[node_id] = np.array(before[node_id])

        # For each node, find the coefficients to map the results.
        model = lambda x, *theta: np.dot(theta, x)

        mapping_coefficients = {}

        for node_id in node_ids_to_map:

            xn, yn = (x[node_id], y[node_id])
            theta, cov = op.curve_fit(model, xn, yn, p0=np.zeros(xn.shape[0]))
            
            mapping_coefficients[node_id] = (theta, cov)

            # Calculate post-mapped values.
            after[node_id] = before[node_id] - model(xn, *theta)

            # A sanity check.
            reference_values = before[node_id] - yn
            before_stddev = np.std(before[node_id] - reference_values)
            after_stddev = np.std(after[node_id] - reference_values)

            assert after_stddev <= before_stddev, \
                "You're doing *worse*?! What?!"

        # Fill in the results for the reference node.
        C = mapping_coefficients.values()[0][0].size
        filler_theta = np.hstack([[1], np.zeros(C - 1)])
        filler_cov = np.zeros((C, C))

        mapping_coefficients[reference_node_id] = (filler_theta, filler_cov)

        return mapping_coefficients


    def _get_mapping_coefficients(self, vector_terms, reference_node_id, **kwargs):

        if 0 > reference_node_id or reference_node_id is None:
            return self._get_mapping_coefficients_from_calibrators(
                vector_terms, **kwargs)

        else:
            return self._get_mapping_coefficients_from_reference_node(
                vector_terms, reference_node_id, **kwargs)


    def _prepare_mapped_data(self, vector_terms, reference_node_id=None,
        default_sigma_calibrator=10e6, minimum_node_estimates=1,
        sql_constraint_for_mapping_query=None,
        sql_constraint_for_data_query=None, **kwargs):
        """
        Prepare data for the noise model so that it can be supplied to Stan.

        This will retrieve the relevant data from the database, and map it on to
        a common reference node (supplied by `reference_node_id`) or to the
        benchmark stars.

        :param vector_terms:
            A list of two-length tuples containing: (1) the parameter(s) to use
            as a single vector term, and (2) the power for that term.

            For example:

            >>> vector_terms=[
                (("TEFF", ), 2),
                (("TEFF", "LOGG"), 1),
                (("LOGG", ), 3.5)
            ]

            Would create a vector term with a scaling offset (always present) 
            and terms of `\theta_1 * TEFF^2`, `\theta_2 * TEFF * LOGG`, and
            `\theta_3 * LOGG^3`.

        :param reference_node_id: [optional]
            The identifier of the reference node to use. If `None` is supplied,
            then the data from all nodes will be mapped to the benchmark star
            results.
    
        :param default_sigma_calibrator: [optional]
            The default uncertainty in calibrator values to use if none exists.

        :param minimum_node_estimates: [optional]
            The minimum number of nodes required for the results to be used in
            the noise model.

        :param sql_constraint_for_mapping_query: [optional]
            Specify an optional SQL constraint to include when retrieving any
            results to be mapped onto a reference node, or the benchmark values.

        :param sql_constraint_for_data_query: [optional]
            Specify an optional SQL constraint to include when retrieving the
            node results that are supplied to the Stan model.
        """

        # Get the mapping coefficients.
        coeffs = self._get_mapping_coefficients(vector_terms, reference_node_id,
            sql_constraint=sql_constraint_for_mapping_query)

        # Construct the query.
        # TECHDEBT: Restricting this to only be able to homogenise teff, logg, feh, xi since they are hard-coded instead of {parameter}

        sql_query = \
            """ WITH    n AS (
                            SELECT id FROM nodes WHERE wg = {wg}), 
                        s AS (
                            SELECT cname, ges_type, ges_fld
                            FROM spectra
                            WHERE ges_type LIKE 'GE_SD_B%') 
                SELECT DISTINCT ON (r.filename, r.node_id) 
                    s.cname, s.ges_type, s.ges_fld,
                    r.filename, r.node_id, r.snr,
                    r.teff, r.e_teff,
                    r.logg, r.e_logg,
                    r.feh, r.e_feh,
                    r.xi, r.e_xi
                FROM    s, n, results as r 
                WHERE r.cname = s.cname 
                  AND r.node_id = n.id 
                  AND r.passed_quality_control = true {sql_constraint}
                """.format(
                    wg=self._wg, parameter=self._parameter,
                    sql_constraint=""   if sql_constraint_for_data_query is None \
                                        else " AND {}".format(sql_constraint_for_data_query))

        # Get the data and apply the mapping coefficients.
        data = self._apply_mapping_coefficients(sql_query, coeffs, vector_terms)

        # Calibrator parameter names
        calibrator_name, calibrator_e_name = _guess_parameter_name(
            self._calibrators, self._parameter)

        finite_calibrator = np.isfinite(self._calibrators[calibrator_name])
        if not np.all(finite_calibrator):
            logger.warn("Not all calibrator values of {} are finite! ({}/{})"\
                .format(self._parameter, sum(finite_calibrator), len(finite_calibrator)))
       
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

                estimates[c, n, v] = data[self._parameter][k]

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

            # Only keep visits that have at least the number of minimum node
            # measurements
            for c in range(C):
                mask = np.sum(np.isfinite(estimates[c]), axis=0) \
                    >= minimum_node_estimates
                
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

        alpha_bounds = dict(teff=(100, 1000), logg=(0.01, 1.0), feh=(0.01, 1.0))
        data_dict.update(
            dict(zip(
                ("lower_alpha_sq", "upper_alpha_sq"), 
                np.array(alpha_bounds[self._parameter])**2)))

        bounds = dict(teff=(3000, 8000), logg=(0, 5), feh=(-3.5, 0.5))
        data_dict.update(
            dict(zip(("lower_bound", "upper_bound"), bounds[self._parameter])))

        # Create additional metadata
        node_names = self._database.retrieve_table("SELECT * FROM nodes")
        metadata = {
            "mapping_coefficients": coeffs,
            "vector_terms": vector_terms,
            "reference_node_id": reference_node_id,
            "calibrators": calibrators,
            "node_ids": unique_estimators,
            "node_names": \
                [node_names["name"][node_names["id"] == node_id][0].strip() \
                    for node_id in unique_estimators]
        }
        result = (data_dict, metadata)
        self._data, self._metadata = result
        return result




    def _apply_mapping_coefficients(self, sql_query, coefficients, vector_terms):
        """
        Apply mapping coefficients to data returned from a SQL query. The
        SQL query must contain all the required parameters outlined in
        `vector_terms`, as well as the `node_id`.
        """

        data = self._database.retrieve_table(sql_query)
        assert data is not None

        # What columns do we require data to have?
        assert "node_id" in data.keys(), \
            "Returned data must contain node_id column"
        for terms, power in vector_terms:
            for term in terms:
                assert term in data.keys(), \
                "Returned data must contain {} column".format(term)

        for node_id, (theta, cov) in coefficients.items():

            match = (data["node_id"] == node_id)

            # Set up the matrix entries.
            x = [np.ones(match.sum())]
            for terms, power in vector_terms:
                value = np.ones(match.sum())
                for term in terms:
                    value *= data[term][match]
                value **= power
                x.append(value)

            # Apply the correction.
            data[self._parameter][match] -= np.dot(theta.T, np.array(x))

        return data



    def _prepare_data(self, default_sigma_calibrator, minimum_node_estimates=1, 
        sql_constraint=None):
        """
        Prepare dthe data for the model so that it can be supplied to Stan.
        
        :param default_sigma_calibrator:
            The default uncertainty in calibrator values to use if none exists.

        :param minimum_node_estimates:
            The minimum number of nodes required for the results to be used in
            the noise model.

        :param sql_constraint: [optional]
            Specify an optional SQL constraint to include when retrieving the
            node results.
        """

        raise DeprecatedError

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
                    wg=self._wg, parameter=self._parameter,
                    sql_constraint=""   if sql_constraint is None \
                                        else " AND {}".format(sql_constraint)))
        assert data is not None, "No calibrator data from WG {}".format(wg)

        # Calibrator parameter names
        calibrator_name, calibrator_e_name = _guess_parameter_name(
            self._calibrators, self._parameter)

        finite_calibrator = np.isfinite(self._calibrators[calibrator_name])
        if not np.all(finite_calibrator):
            logger.warn("Not all calibrator values of {} are finite! ({}/{})"\
                .format(self._parameter, sum(finite_calibrator), len(finite_calibrator)))
       
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

                estimates[c, n, v] = data[self._parameter][k]

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

            # Only keep visits that have at least the number of minimum node
            # measurements
            for c in range(C):
                mask = np.sum(np.isfinite(estimates[c]), axis=0) \
                    >= minimum_node_estimates
                
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
            dict(zip(("lower_alpha_sq", "upper_alpha_sq"), np.array(alpha_bounds[self._parameter])**2)))

        bounds = dict(teff=(3000, 8000), logg=(0, 5), feh=(-3.5, 0.5))
        data_dict.update(
            dict(zip(("lower_bound", "upper_bound"), bounds[self._parameter])))

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