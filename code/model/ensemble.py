
"""
Classes to deal with homogenisation models written in Stan.
"""

import cPickle as pickle
import logging
import numpy as np
import os
import pystan as stan
from astropy.table import Table
from astropy.io import fits # HACK REMOVE TODO


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


class BaseEnsembleModel(object):

    def __init__(self, database, wg, calibrators, recompile=False, 
        overwrite=True, **kwargs):
        self._database = database
        self._wg = wg
        self._calibrators = calibrators

        self._model = self._load_model(recompile, overwrite, **kwargs)
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

        return kwds


    def optimize(self, data, recompile=False, overwrite=False, **kwargs):
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


    def sample(self, data, chains=4, iter=2000, warmup=None, recompile=False, 
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

        return self._model.sampling(**kwds)


class NewEnsembleModel(BaseEnsembleModel):

    _MODEL_PATH = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "new-ensemble-model.stan") # new-ensemble-model.stan

    def _prepare_data(self, parameter, default_sigma_calibrator=1e3,
        fill_function=np.mean, fill_variance=1e10):
        """
        Prepare the data for the model so that it can be supplied to STAN.

        :param parameter:
            The name of the model parameter (e.g., teff) that will be used in
            this single parameter ensemble model.

        :param fill_function: [optional]
            A function that takes as input finite values and returns as output a
            value that will be used to fill arrays that would otherwise be non-
            finite entries. As long as the `fill_variance` is suitably large, it
            does not matter what the fill value is.

        :param fill_variance: [optional]
            The variance to add for matrix entries that have been filled with
            the `fill_function` value.
        """


        parameter = str(parameter).lower()
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
                """.format(wg=self._wg, parameter=parameter))
        assert data is not None, "No calibrator data from WG {}".format(wg)

        # Calibrator parameter name
        calibrator_name, calibrator_e_name = _guess_parameter_name(
            self._calibrators, parameter)

        finite_calibrator = np.isfinite(self._calibrators[calibrator_name])
        if not np.all(finite_calibrator):
            logger.warn("Not all calibrator values of {} are finite! ({}/{})"\
                .format(parameter, sum(finite_calibrator), len(finite_calibrator)))

        calibrators = self._calibrators[finite_calibrator]

        N_calibrators = sum(finite_calibrator)
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
            spectrum_snr[i] = data["snr"][si]
            assert np.std(data["snr"][si:ei]) < 1

            for k in range(si, ei):
                l = np.where(data["node_id"][k] == unique_estimators)[0][0]
                estimates[i, l] = data[parameter][k]
        
        # Fill in the NaN estimates with the mean of that entry, and create an
        # additive variance array.
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
        keep = np.all(np.isfinite(estimates), axis=1)
        
        # Ensure spectrum_snr values are valid.
        spectrum_snr[spectrum_snr < 1] = 1
        data_dict.update({
            "ivar_spectrum": (1.0/spectrum_snr[keep])**2,
            "N_calibrator_visits": sum(keep),
            "var_additive": var_additive[keep],
            "estimates": estimates[keep],
            # Make the calibrator index 1-indexed for STAN
            "calibrator_index": 1 + calibrator_index[keep].astype(int),
        })

        # Update with additional metadata
        node_names = self._database.retrieve_table("SELECT * FROM nodes")
        raise a
        data_dict.update({
            "__N_pairwise_estimators": np.math.factorial(N_estimators) \
                / (2*np.math.factorial(max(1, N_estimators - 2))),
            "__calibrator_ges_fld": unique_calibrators,
            "__node_ids": unique_estimators,
        #    "__node_names":
        })

        return data_dict

        
















class SingleParameterEnsembleModel(BaseEnsembleModel):

    """
    An ensemble model for inferring the accuracy and precision of different
    nodes (estimators) of a single parameter.

    :param database:
        A database to connect to.

    :param wg:
        The working group number (e.g., 11).

    :param calibrators:
        A table of stars with previously determined parameters, or additional
        information separate from spectroscopy (e.g., asteroseismic surface
        gravities, isochrones, etc). This table must include a `GES_FLD`
        column (and eventually should accept `CNAME` columns) to match with
        against the database.

    :param recompile: [optional]
        Re-compile the model if it has already been compiled. Default is False.

    :param overwrite: [optional]
        Overwrite the compiled model path if it already exists. Default is False.
    """

    _MODEL_PATH = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "ensemble-model.stan")

    # MAGIC NUMBER: default_calibrator_sigma, fill_variance
    def _prepare_data(self, parameter, default_calibrator_sigma=1e3, 
        fill_function=np.mean, fill_variance=1e10):
        """
        Prepare the data for the model so that it can be supplied to STAN.

        :param parameters:
            The name of the model parameter (e.g., teff) that will be used in
            this single parameter ensemble model.

        :param fill_function: [optional]
            A function that takes as input finite values and returns as output a
            value that will be used to fill arrays that would otherwise be non-
            finite entries. As long as the `fill_variance` is suitably large, it
            does not matter what the fill value is.

        :param fill_variance: [optional]
            The variance to add for matrix entries that have been filled with
            the `fill_function` value.
        """

        parameter = str(parameter).lower()
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
                    r.filename, r.node_id, s.cname, s.ges_type, s.ges_fld,
                    r.{parameter}, r.e_{parameter}
                FROM    s, n, results as r 
                WHERE r.cname = s.cname 
                  AND r.node_id = n.id 
                """.format(wg=self._wg, parameter=parameter))
        assert data is not None, "No calibrator data from WG {}".format(wg)

        # Calibrator parameter name
        calibrator_name, calibrator_e_name = _guess_parameter_name(
            self._calibrators, parameter)

        finite_calibrator = np.isfinite(self._calibrators[calibrator_name])
        if not np.all(finite_calibrator):
            logger.warn("Not all calibrator values of {} are finite! ({}/{})"\
                .format(parameter, sum(finite_calibrator), len(finite_calibrator)))

        calibrators = self._calibrators[finite_calibrator]

        N_calibrators = sum(finite_calibrator)
        N_estimators = len(set(data["node_id"]))
        data_dict = {
            "N_calibrators": N_calibrators,
            "N_estimators": N_estimators,
            "calibrator_mu": np.array(calibrators[calibrator_name]),
            "calibrator_sigma": np.array(calibrators[calibrator_e_name]),
        }

        # Check the data so far.
        if 1 > N_calibrators:
            raise ValueError("no calibrators")

        if 1 > N_estimators:
            raise ValueError("no estimators")

        if not np.all(np.isfinite(data_dict["calibrator_sigma"])):
            logger.warn(
                "Not all finite calibrator entries of {} have finite errors! " \
                "Setting default calibrator_sigma as {}".format(
                    parameter, default_calibrator_sigma))
            finite = np.isfinite(data_dict["calibrator_sigma"])
            data_dict["calibrator_sigma"][~finite] = default_calibrator_sigma

        # OK now update the data dictionary with the spectroscopic measurements.
        # Need to group by node id and CNAME.
        
        mean_estimate = np.nan * np.ones((N_calibrators, N_estimators))
        N_estimates_per_estimator = np.zeros(mean_estimate.shape, dtype=int)
        # The unique calibrators needs to be determined by the calibrator
        # identifier, that way we won't be including data from the database 
        # for which we don't have a calibrated value.
        # NOTE: Don't sort the unique_calibrators, because it must align with
        #       the calibrator_mu and calibrator_sigma above.
        unique_calibrators = np.array(map(str.strip, calibrators["GES_FLD"]))
        unique_estimators = np.sort(np.unique(data["node_id"]))
        assert len(set(calibrators["GES_FLD"])) == N_calibrators
        assert len(unique_estimators) == N_estimators

        # Group data by the GES_FLD and node_id
        skipped = []
        data = data.group_by(["ges_fld", "node_id"])
        for i, si in enumerate(data.groups.indices[:-1]):
            ei = data.groups.indices[i + 1]

            ges_fld = data["ges_fld"][si].strip()

            try:
                j = np.where(ges_fld == unique_calibrators)[0][0]

            except IndexError:
                # There is a calibrator in the database that we're not using here
                if ges_fld not in skipped:
                    logger.warn("Skipping calibrator {}".format(ges_fld))
                    skipped.append(ges_fld)
                continue

            k = np.where(data["node_id"][si] == unique_estimators)[0][0]

            # Calculate the average parameter value of the group.
            estimates = data[parameter][si:ei]
            mean_estimate[j, k] = np.nanmean(estimates)
            
            # And the number of estimates that went into this mean.
            N_estimates_per_estimator[j, k] = np.sum(np.isfinite(estimates))
            
        # Check that things make sense: we should only have NaNs if there were
        # zero estimates for a given value.
        assert np.any(np.isfinite(mean_estimate))
        assert np.all(
            ~np.isfinite(mean_estimate) == (N_estimates_per_estimator == 0))

        # Create additive variance array and fill the mean estimates with some
        # fill value.
        no_estimates = ~np.isfinite(mean_estimate)
        assert not np.all(no_estimates)
        mean_estimate[no_estimates] = fill_function(mean_estimate[~no_estimates])

        additive_var = np.zeros_like(mean_estimate)
        additive_var[no_estimates] += fill_variance

        # We also need to update the number of estimates array so that this is
        # considered to be a single (bad) estimate, otherwise we will divide by
        # zero.
        N_estimates_per_estimator[no_estimates] = 1

        # Update the data dict.
        data_dict.update({
            "mean_estimate": mean_estimate,
            "N_estimates_per_estimator": N_estimates_per_estimator,
            "additive_var": additive_var,
            "__N_pairwise_estimators": np.math.factorial(N_estimators)/(2*np.math.factorial(max(1, N_estimators - 2))),
            
        })

        assert np.all([np.all(np.isfinite(v)) for v in data_dict.values()])

        return data_dict


class SingleParameterEnsembleModelWithCorrelations(SingleParameterEnsembleModel):

    _MODEL_PATH = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "single-parameter-ensemble-model-with-correlations.stan")


class MultipleParameterEnsembleModel(BaseEnsembleModel):

    """
    An ensemble model for inferring the accuracy and precision of different
    nodes (estimators) of multiple parameters, simultaneously.

    :param database:
        A database to connect to.

    :param wg:
        The working group number (e.g., 11).

    :param calibrators:
        A table of stars with previously determined parameters, or additional
        information separate from spectroscopy (e.g., asteroseismic surface
        gravities, isochrones, etc). This table must include a `GES_FLD`
        column (and eventually should accept `CNAME` columns) to match with
        against the database.

    :param recompile: [optional]
        Re-compile the model if it has already been compiled. Default is False.

    :param overwrite: [optional]
        Overwrite the compiled model path if it already exists. Default is False.
    """

    _MODEL_PATH = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "multiple-parameter-ensemble-model.stan")

    # MAGIC NUMBER: default_calibrator_sigma, fill_variance
    def _prepare_data(self, parameters, default_calibrator_sigmas, 
        fill_function=np.mean, fill_variance=1e10):
        """
        Prepare the data for the model so that it can be supplied to STAN.

        :param parameters:
            A list-like object containing the name of the model parameters that
            will be used in this ensemble model.

        :param default_calibrator_sigmas:
            A list-like object containing the default uncertainties to adopt in
            each parameter for calibration objects without an uncertainty.

        :param fill_function: [optional]
            A function that takes as input finite values and returns as output a
            value that will be used to fill arrays that would otherwise be non-
            finite entries. As long as the `fill_variance` is suitably large, it
            does not matter what the fill value is.

        :param fill_variance: [optional]
            The variance to add for matrix entries that have been filled with
            the `fill_function` value.
        """

        for parameter in parameters:
            parameter = str(parameter).lower()
            if parameter not in ("teff", "logg", "feh", "xi", "mh"):
                raise AreYouSureYouKnowWhatYoureDoing

        N_parameters = len(parameters)

        # Get the data from the database for this WG.
        # TODO: Need a better way to identify calibrators. Right now we do it
        #       just on GES_TYPE, but in future we may want to do this directly
        #       on to CNAMEs in the calibrator list.
        parameter_sql = ", ".join(["r.{0}, r.e_{0}".format(p) for p in parameters])
        data = self._database.retrieve_table(
            """ WITH    n AS (
                            SELECT id FROM nodes WHERE wg = {wg}), 
                        s AS (
                            SELECT cname, ges_type, ges_fld
                            FROM spectra
                            WHERE ges_type LIKE 'GE_SD_B%') 
                SELECT DISTINCT ON (r.filename, r.node_id) 
                    r.filename, r.node_id, s.cname, s.ges_type, s.ges_fld,
                    {parameter_sql}
                FROM    s, n, results as r 
                WHERE r.cname = s.cname 
                  AND r.node_id = n.id 
                """.format(wg=self._wg, parameter_sql=parameter_sql))
        assert data is not None, "No calibrator data from WG {}".format(wg)

        # Construct arrays of calibrator_mu and calibrator_sigma, then check for
        # finite entries.

        # Calibrator parameter name
        calibrator_names, calibrator_e_names \
            = _guess_parameter_names(self._calibrators, parameters)

        calibrator_mu \
            = np.vstack([self._calibrators[n] for n in calibrator_names]).T
        calibrator_sigma \
            = np.vstack([self._calibrators[n] for n in calibrator_e_names]).T

        finite_calibrator_mu = np.all(np.isfinite(calibrator_mu), axis=0)
        for index in np.where(~finite_calibrator_mu)[0]:
            logger.warn("Not all calibrator values of {} ({}) are finite!"\
                .format(parameters[index], calibrator_names[index]))

        finite_calibrator_sigma = np.all(np.isfinite(calibrator_sigma), axis=0)
        for index in np.where(~finite_calibrator_sigma)[0]:

            not_finite = ~np.isfinite(calibrator_sigma[:, index])
            logger.warn(
                "Not all calibrator values of {} ({}) are finite: "\
                "{}/{} not finite: {}".format(
                    parameters[index], calibrator_e_names[index],
                    sum(not_finite), len(not_finite),
                    ", ".join([ges_fld.strip() \
                        for ges_fld in self._calibrators["GES_FLD"][~not_finite]])))

            calibrator_sigma[not_finite, index] \
                = default_calibrator_sigmas[index]


        finite = np.all(np.isfinite(calibrator_mu * calibrator_sigma), axis=1)
        calibrators = self._calibrators[finite]
        calibrator_mu = calibrator_mu[finite]
        calibrator_sigma = calibrator_sigma[finite]

        N_calibrators = sum(finite)
        N_estimators = len(set(data["node_id"]))
        data_dict = {
            "N_calibrators": N_calibrators,
            "N_parameters": N_parameters,
            "N_estimators": N_estimators,
            "calibrator_mu": calibrator_mu.T.flatten(),
            "calibrator_sigma": calibrator_sigma.T.flatten()
        }

        # Check the data so far.
        if 1 > N_calibrators:
            raise ValueError("no calibrators")

        if 1 > N_estimators:
            raise ValueError("no estimators")

        # OK now update the data dictionary with the spectroscopic measurements.
        # Need to group by node id and CNAME.
        mean_estimate = np.nan * np.ones((N_calibrators, N_estimators * N_parameters))
        N_estimates_per_estimator = np.zeros(mean_estimate.shape, dtype=int)
        # The unique calibrators needs to be determined by the calibrator
        # identifier, that way we won't be including data from the database 
        # for which we don't have a calibrated value.
        # NOTE: Don't sort the unique_calibrators, because it must align with
        #       the calibrator_mu and calibrator_sigma above.
        unique_calibrators = np.array(map(str.strip, calibrators["GES_FLD"]))
        unique_estimators = np.sort(np.unique(data["node_id"]))
        assert len(set(calibrators["GES_FLD"])) == N_calibrators
        assert len(unique_estimators) == N_estimators

        # Group data by the GES_FLD and node_id
        skipped = []
        data = data.group_by(["ges_fld", "node_id"])
        for idx, si in enumerate(data.groups.indices[:-1]):
            ei = data.groups.indices[idx + 1]

            ges_fld = data["ges_fld"][si].strip()

            try:
                i = np.where(ges_fld == unique_calibrators)[0][0]

            except IndexError:
                # There is a calibrator in the database that we're not using here
                if ges_fld not in skipped:
                    logger.warn("Skipping calibrator {}".format(ges_fld))
                    skipped.append(ges_fld)
                continue

            j = np.where(data["node_id"][si] == unique_estimators)[0][0]

            # Calculate the average parameter value of the group.
            for k, parameter in enumerate(parameters):
                estimates = data[parameter][si:ei]
                mean_estimate[i, j + k*N_estimators] = np.nanmean(estimates)
            
                # And the number of estimates that went into this mean.
                N_estimates_per_estimator[i, j + k*N_estimators] = np.sum(np.isfinite(estimates))
            
        # Check that things make sense: we should only have NaNs if there were
        # zero estimates for a given value.
        assert np.any(np.isfinite(mean_estimate))
        assert np.all(
            ~np.isfinite(mean_estimate) == (N_estimates_per_estimator == 0))

        # Create additive variance array and fill the mean estimates with some
        # fill value.
        additive_var = np.zeros_like(mean_estimate)
        for k in range(len(parameters)):
            sj, ej = (N_estimators * k, N_estimators * (k + 1))

            no_estimates = ~np.isfinite(mean_estimate[:, sj:ej])

            assert not np.all(no_estimates)
            # Run the fill function on all values for that particular parameter
            mean_estimate[:, sj:ej][no_estimates] \
                = fill_function(mean_estimate[:, sj:ej][~no_estimates])
            additive_var[:, sj:ej][no_estimates] += fill_variance

        # We also need to update the number of estimates array so that this is
        # considered to be a single (bad) estimate, otherwise we will divide by
        # zero.
        N_estimates_per_estimator[N_estimates_per_estimator == 0] = 1

        # Update the data dict.
        data_dict.update({
            "mean_estimate": mean_estimate,
            "N_estimates_per_estimator": N_estimates_per_estimator,
            "additive_var": additive_var,
            "additive_sigma": additive_var**0.5,
            # Calculate the number of pairwise coefficients?
            "__N_pairwise_estimators": N_parameters * np.math.factorial(N_estimators)/(2*np.math.factorial(max(1, N_estimators - 2))),
            "__N_pairwise_parameters": np.math.factorial(N_parameters)/(2*np.math.factorial(max(1, N_parameters - 2)))
        })
        data_dict["sqrt_N_estimates_per_estimator"] = N_estimates_per_estimator**0.5

        data_dict["__node_ids"] = unique_estimators
        #data_dict["__calibrator_ges_fld"] = unique_calibrators

        return data_dict


# Usage
"""
ok = benchmarks["TEFF"] < 8000
homogenisation_model = EnsembleModel(database, wg, benchmarks[ok])

homogenisation_model.optimize()
homogenisation_model.sample()
homogenisation_model.plot_trace()
homogenisation_model.plot_corner()

homogenisation_model.homogenize(cname=XX, exclude=??)
"""


