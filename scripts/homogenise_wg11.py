
"""
Homogenisation models.
"""

import yaml
import logging
import numpy as np

from code import GESDatabase
from code.model.ensemble import  SingleParameterEnsembleModel

from astropy.table import Table

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

wg = 11
parameters = ("teff", "logg")
sample_kwds = dict(chains=4, iter=2000)

scales = dict(teff=250, logg=0.25)

models = {}
for parameter in parameters:
    model = SingleParameterEnsembleModel(database, wg, parameter, benchmarks)

    scale = scales[parameter]
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

    model.write(
        "homogenisation-wg{}-{}.model".format(wg, parameter), overwrite=True)

    # Homogenise this parameter for all stars.
    model.homogenise_all_stars(update_database=True)

    models[parameter] = model