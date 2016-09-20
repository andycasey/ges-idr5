
"""
Sandbox for ensemble model
"""

import yaml
import logging
import os
import matplotlib.pyplot as plt
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





from code.model.ensemble import SingleParameterEnsembleModel

benchmarks = Table.read("fits-templates/benchmarks/GES_iDR5_FGKMCoolWarm_Benchmarks_AcceptedParams_01082016.fits")

model = SingleParameterEnsembleModel(database, 11, benchmarks[benchmarks["TEFF"] < 8000])


data_dict = model._prepare_data("teff", default_calibrator_sigma=150)

op_params = model.optimize(data_dict, init={
    "var_intrinsic": 100**2,
    "var_node_rand": np.ones(data_dict["N_estimators"]) * 100**2,
    "var_node_sys": np.ones(data_dict["N_estimators"]) * 100**2,
    "truths": np.array(data_dict["calibrator_mu"])
    })

# Drop the optimized covariance matrix.
del op_params["covariance"]

samples = model.sample(data_dict, init=op_params)

fig = samples.plot(pars=("var_intrinsic", "var_node_sys", "var_node_rand"))