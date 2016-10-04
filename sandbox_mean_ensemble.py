
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





from code.model.ensemble import  MeanEnsembleModel

benchmarks = Table.read("fits-templates/benchmarks/GES_iDR5_FGKMCoolWarm_Benchmarks_AcceptedParams_01082016.fits")
benchmarks["E_FEH"] = 0.1


model = MeanEnsembleModel(database, 11, "teff", benchmarks[benchmarks["TEFF"] < 8000])

data, metadata = model._prepare_data()

init = {
    "truths": data["mu_calibrator"],
    "var_sys_estimator": 100**2 * np.ones(data["N_nodes"]),
    "alpha_sq": np.ones(data["N_nodes"]),
    "rho_estimators": np.zeros(metadata["N_pairwise_nodes"]),
    "bias": np.zeros(data["N_nodes"])
}

op_params = model.optimize(data, init=init)

fit = model.sample(data, init=op_params, iter=10)
raise a

model.homog






raise a


