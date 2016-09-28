
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

model = NewEnsembleModel(database, 11, benchmarks[benchmarks["TEFF"] < 8000])

data = model._prepare_data("teff", )

op_params = model.optimize(data, init={
    "truths": data["mu_calibrator"],
    "var_intrinsic": 100**2,
    "var_sys_estimator": 100**2 * np.ones(7),
    "alpha_sq": 1000 * np.ones(7),
    })

samples = model.sample(data, init=op_params)








raise a