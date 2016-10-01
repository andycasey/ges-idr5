
"""
Homogenisation models.
"""

import yaml
import logging
import numpy as np
import os
from astropy.table import Table
from collections import OrderedDict


from code import GESDatabase
from code.model.ensemble import  SingleParameterEnsembleModel

# Initialize logging.
logger = logging.getLogger("ges")

# Create a database object.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)

# Load the "benchmarks"
# Only use "benchmarks" with TEFF < 8000 K
benchmarks = Table.read("fits-templates/benchmarks/GES_iDR5_FGKM_Benchmarks_ARC_29092016.fits")
benchmarks = benchmarks[benchmarks["TEFF"] < 8000]
benchmarks["E_FEH"] = 0.15

model_paths = "homogenisation-wg{wg}-{parameter}.model"

wgs = (11, 12, 10)
parameter_scales = OrderedDict([
    ("teff", 250),
    ("logg", 0.25),
    ("feh", 0.25)
])

sample_kwds = dict(chains=4, iter=2000)

models = {}
for wg in wgs:
    models[wg] = {}
    for parameter, scale in parameter_scales.items():

        model_path = model_paths.format(wg=wg, parameter=parameter)

        if os.path.exists(model_path):
            model = SingleParameterEnsembleModel.read(model_path, database)

        else:
            model = SingleParameterEnsembleModel(database, wg, parameter, benchmarks)
            data, metadata = model._prepare_data(default_sigma_calibrator=scale)

            init = {
                "truths": data["mu_calibrator"],
                "var_intrinsic": scale**2,
                "var_sys_estimator": scale**2 * np.ones(data["N_estimators"]),
                "alpha_sq": 1000 * np.ones(data["N_estimators"]),
                "rho_estimators": np.zeros(metadata["N_pairwise_estimators"]),
                "c0_estimators": np.zeros(data["N_estimators"])
            }

            op_params = model.optimize(data, init=init)

            fit = model.sample(data, init=op_params, **sample_kwds)

            model.write(model_path, overwrite=True)

        # Homogenise this parameter for all stars.
        #model.homogenise_all_stars(update_database=True)

        # Keep the model.
        models[wg][parameter] = model

# TODO:
# - VROT
# - VEL
# - FLAGS
# - XI

# Some notes:
# - If a node provides `mh` or `feh`, we can treat them as the same (see `scripts/setup_db.py`)
