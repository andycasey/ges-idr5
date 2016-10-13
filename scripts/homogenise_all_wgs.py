
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
from code.model.ensemble import  SingleParameterEnsembleModel, EnsembleModelWithSysVarianceModel

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
benchmarks["E_FEH"] = 0.10

model_paths = "homogenisation-wg{wg}-{parameter}.model"

wgs = (11, )
parameter_scales = OrderedDict([
    ("teff", 250),
    ("logg", 0.25),
    ("feh", 0.25)
])

sample_kwds = dict(chains=4, iter=2000)

finite = np.isfinite(benchmarks["TEFF"] * benchmarks["LOGG"] * benchmarks["FEH"])
benchmarks = benchmarks[finite]

models = {}
for wg in wgs:
    models[wg] = {}
    for parameter, scale in parameter_scales.items():

        model_path = model_paths.format(wg=wg, parameter=parameter)

        if os.path.exists(model_path):
            model = SingleParameterEnsembleModel.read(model_path, database)

        else:
            model = EnsembleModelWithSysVarianceModel(database, wg, parameter, benchmarks)
            data, metadata = model._prepare_data(default_sigma_calibrator=scale)

            init = {
                # -1 offset to account for Stan's 1-indexing.
                "truths": data["mu_calibrator"],
                #"var_sys_estimator": (scale/5.)**2 * np.ones(data["N_estimators"]),
                "alpha": np.mean([data["lower_alpha"], data["upper_alpha"]]) \
                    * np.ones(data["N_estimators"]),
                "rho_estimators": np.zeros(metadata["N_pairwise_estimators"]),
                "c0_estimators": np.zeros(data["N_estimators"]),

            }

            init["vs_log_c"] = np.log(np.ones(data["N_estimators"]) * scale)
            init["vs_a"] = np.ones(data["N_estimators"])
            init["vs_b"] = 2 * np.ones(data["N_estimators"])

            raise a
            op_params = model.optimize(data, init=init, iter=100000)

            raise a
            fit = model.sample(data, init=op_params, **sample_kwds)

            model.write(model_path, overwrite=True)

        # Keep the model.
        models[wg][parameter] = model


        # Homogenise this parameter for all stars.
        model.homogenise_all_stars(update_database=True)


# TODO:
# - VROT
# - VEL
# - FLAGS
# - XI

# Some notes:
# - If a node provides `mh` or `feh`, we can treat them as the same (see `scripts/setup_db.py`)
