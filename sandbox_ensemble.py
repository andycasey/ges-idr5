
"""
Sandbox for ensemble model
"""

import yaml
import logging
import os
import matplotlib.pyplot as plt
import numpy as np
from glob import glob

from code import (GESDatabase, plot, summary)
from astropy.table import Table

# Initialize logging.
logger = logging.getLogger("ges.idr5.qc")
logger.setLevel(logging.DEBUG)

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    "%(asctime)s [%(levelname)-8s] %(message)s"))
logger.addHandler(handler)


# Create a database object.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)





from code.model.ensemble import (NewEnsembleModel, SingleParameterEnsembleModel, MultipleParameterEnsembleModel, SingleParameterEnsembleModelWithCorrelations)

benchmarks = Table.read("fits-templates/benchmarks/GES_iDR5_FGKMCoolWarm_Benchmarks_AcceptedParams_01082016.fits")



"""
model = SingleParameterEnsembleModel(database, 11, benchmarks[benchmarks["TEFF"] < 8000])

data = model._prepare_data("teff", default_calibrator_sigma=150)

op_params = model.optimize(data, init={
    "intrinsic_var": 100**2,
    "estimator_sys_var": np.ones(data["N_estimators"]) * 100**2,
    "estimator_rand_var": np.ones(data["N_estimators"]) * 100**2,
    "truths": np.array(data["calibrator_mu"]),
    "rho_parameters": np.zeros(10)
    }, overwrite=True)

# Drop the optimized covariance matrix.
del op_params["covariance"]

samples = model.sample(data, init=op_params, overwrite=True)

fig = samples.plot(
    pars=("intrinsic_sigma", "estimator_rand_sigma", "estimator_sys_sigma"))



"""


benchmarks["E_FEH"] = 0.1


"""
model = MultipleParameterEnsembleModel(database, 11, 
    benchmarks[(8000 > benchmarks["TEFF"]) * (benchmarks["TEFF"] > 4000)])

data = model._prepare_data(("teff", "logg"), (150, 0.1, ))

intrinsic_sigma = np.array([150, 0.1])
K = data["N_parameters"] * data["N_estimators"]
init = {
    #"intrinsic_var": intrinsic_var,
    "estimator_sys_sigma": np.tile(intrinsic_sigma, data["N_estimators"]).reshape(data["N_estimators"], -1),
    "estimator_rand_sigma": np.tile(intrinsic_sigma, data["N_estimators"]).reshape(data["N_estimators"], -1),
    "truths": data["calibrator_mu"].reshape(-1, data["N_calibrators"]).T,
    "rho_estimators": np.zeros(data["__N_pairwise_estimators"]),
    "rho_parameters": np.zeros(data["__N_pairwise_parameters"])
    #"Sigma": 100 * np.array([np.eye(data["N_parameters"] * data["N_estimators"]) for i in range(data["N_calibrators"])])
    }
#op_params = model.optimize(data, iter=10000, init=init)


bar = model.sample(data, init=lambda x=0: init, iter=10000, chains=1, validate=False)
"""


"""
model = SingleParameterEnsembleModelWithCorrelations(database, 11, benchmarks[benchmarks["TEFF"] < 8000])

data = model._prepare_data("teff", default_calibrator_sigma=150)
data["additive_sigma"] = data["additive_var"]**0.5

intrinsic_sigma = 150.0
init = {
    "estimator_sys_sigma": np.tile(intrinsic_sigma, data["N_estimators"]),
    "estimator_rand_sigma": np.tile(intrinsic_sigma, data["N_estimators"]),
    "truths": data["calibrator_mu"],
    "rho_estimators": np.zeros(data["__N_pairwise_estimators"]),
    "Sigma": intrinsic_sigma**2 * np.array([np.eye(data["N_estimators"]) for _ in range(data["N_calibrators"])])
    }

op_params = model.optimize(data, init=init)
"""

model = NewEnsembleModel(database, 10, benchmarks[benchmarks["TEFF"] < 8000])

data, metadata = model._prepare_data("teff")

var_intrinsic = 100**2
init = {
    "truths": data["mu_calibrator"],
    "var_intrinsic": var_intrinsic**2,
    "var_sys_estimator": var_intrinsic**2 * np.ones(data["N_estimators"]),
    "alpha_sq": 1000 * np.ones(data["N_estimators"]),
    "rho_estimators": np.zeros(metadata["N_pairwise_estimators"]),
    "c0_estimators": np.zeros(data["N_estimators"])
}

op_params = model.optimize(data, init=init)

fit = model.sample(data, init=op_params)







raise a




def homogenise_survey_measurements(cname, wg, parameter, ensemble_model_samples,
    database):
    """
    Produce an unbiased estimate of an astrophyiscal parameter for a given
    survey object.

    :param cname:
        The CNAME (unique star identifier) of an object.

    :param wg:
        The working group to consider measurements from.

    :param parameter:
        The name of the parameter to estimate.

    :param ensemble_model_samples:
        Samples from the ensemble model to use when estimating the astrophysical
        parameter.
    """

    # Get the data for this object.
    estimates = database.retrieve_table(
        """ SELECT  DISTINCT ON (filename, node_id)
                    cname, node_id, snr, {parameter}
            FROM    results, nodes
            WHERE   nodes.wg = {wg}
              AND   nodes.id = results.node_id
              AND   cname = '{cname}'
              AND   {parameter} <> 'NaN'
        """.format(wg=wg, cname=cname, parameter=parameter))

    assert estimates is not None

    # Extract N samples for all the parameters.

    # For each sample, calculate:
    #   1. The total variance (systematic**2 + (alpha/SNR)**2)
    #   2. The weighted mean from all observations by that nodes.
    #    --> check that this follows 1/sqrt(N)
    #   3. Construct a covariance matrix using the weighted means, uncertainties
    #       and the correlation coefficients
    #   4. Draw from a Gaussian using the weighted means and your new Cov matrix
    #   5. Record the draw.

    pars = [
        "var_intrinsic",
        "var_sys_estimator",
        "alpha_sq",
        "rho_estimators",
        "c0_estimators"
    ]

    
    samples = ensemble_model_samples.extract(pars=pars)

    unique_node_ids = ensemble_model_samples.data["node_ids"]
    K = len(samples["var_intrinsic"])

    estimates = estimates.group_by("node_id")
    
    # 1. Calculate the total variance in each measurement.
    var_total = np.zeros((len(estimates), K))
    for j in range(len(estimates)):

        # Get the node index.
        k = np.where(estimates["node_id"][j] == unique_node_ids)[0][0]

        var_total[j, :] \
            = samples["var_sys_estimator"][:, k] \
                + samples["alpha_sq"][:, k]/estimates["snr"][j]

                
    # 2. Calculate the weighted mean from each node.
    M = len(set(estimates["node_id"]))
    weighted_node_mu = np.zeros((M, K))
    weighted_node_variance = np.zeros((M, K))
    node_ids = np.zeros(M)
    for i, si in enumerate(estimates.groups.indices[:-1]):
        ei = estimates.groups.indices[i + 1]

        mu = (estimates[parameter][si:ei]).reshape(-1, 1) # Biases
        variance = var_total[si:ei]

        weights = 1.0/variance
        normalized_weights = weights/np.sum(weights, axis=0)


        weighted_mu = np.sum(normalized_weights * mu, axis=0)
        weighted_variance = 1.0/np.sum(weights, axis=0)

        weighted_node_mu[i, :] = weighted_mu + samples["c0_estimators"][:, i]
        weighted_node_variance[i, :] = weighted_variance
        node_ids[i] = estimates["node_id"][si]

    posterior = np.nan * np.ones(K)
    for i in range(K):

        Sigma = np.eye(M) * weighted_node_variance[:, i]
        
        a = 0
        for j in range(M):
            for k in range(j + 1, M):
                term = samples["rho_estimators"][i, a] * Sigma[j, j]**0.5 * Sigma[k, k]**0.5
                Sigma[j, k] = term
                Sigma[k, j] = term
                a += 1

        W = np.ones((M, 1))
        Cinv = np.linalg.inv(Sigma)
        var_min = 1.0/np.dot(np.dot(W.T, Cinv), W)
        posterior[i] = var_min * np.dot(np.dot(W.T, Cinv), weighted_node_mu[:, i])
    
    return posterior




sun = homogenise_survey_measurements(
    "ssssssss-sssssss", 11, "logg", samples, database)




# 18 Sco:
homogenise_survey_measurements(
    "16153746-0822162", 11, "logg", samples, database)






homogenise_survey_measurements(
    "14153967+1910567", 11, "logg", samples, database)