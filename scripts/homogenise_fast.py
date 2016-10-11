
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
from code.model.ensemble import EnsembleModel, MedianModel

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

sample_kwds = dict(chains=4, iter=10000, thin=10)

finite = np.isfinite(benchmarks["TEFF"] * benchmarks["LOGG"] * benchmarks["FEH"])
benchmarks = benchmarks[finite]

model = MedianModel(database, 11, "xi")
model.homogenise_all_stars(update_database=True, default_sigma=0.5)

model = MedianModel(database, 11, "alpha_fe")
model.homogenise_all_stars(update_database=True, default_sigma=0.10)


models = {}
for wg in wgs:
    models[wg] = {}
    for parameter, scale in parameter_scales.items():

        model_path = model_paths.format(wg=wg, parameter=parameter)
        if os.path.exists(model_path): 
            model = EnsembleModel.read(model_path, database)

        else:
            model = EnsembleModel(database, wg, parameter, benchmarks)
            data, metadata = model._prepare_data(
                default_sigma_calibrator=scale)

            init = {
                "truths": data["mu_calibrator"],
                "biases": np.zeros(data["N"]),
                "missing_estimates": np.random.uniform(
                    data["lower_bound"], data["upper_bound"], size=data["TM"]),
                "alpha_sq": np.mean([data["lower_alpha_sq"], data["upper_alpha_sq"]]) * np.ones(data["N"]),
                "vs_c": scale**2 * np.ones(data["N"]),
                "vs_a": 1e-2 + np.zeros((data["N"], data["S"])).T,
                "vs_b": 1e-2 + np.ones((data["N"], data["S"])).T,
                "L_corr": np.eye(data["N"])
            }

            print("Number of model parameters: {}".format(
                sum([np.array(v).size for v in init.values()])))

            op_params = model.optimize(data, init=init, iter=100000)

            fit = model.sample(data, init=op_params, **sample_kwds)

            model.write(model_path, 
                overwrite=True, __ignore_model_pars=("Sigma", "full_rank_estimates"))


        model.homogenise_all_stars(update_database=True)


# TODO:
# - VROT
# - VEL
# - FLAGS
# - XI

# Some notes:
# - If a node provides `mh` or `feh`, we can treat them as the same (see `scripts/setup_db.py`)
