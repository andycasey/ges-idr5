
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
from code.model.ensemble import MissingDataModel

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

        model = MissingDataModel(database, wg, parameter, benchmarks)
        data, metadata = model._prepare_data(
            default_sigma_calibrator=scale)
        

        init = {
            "truths": data["mu_calibrator"],
            "biases": np.zeros(data["N"]),
            "missing_estimates": np.mean([data["upper_bound"], data["lower_bound"]]) * np.ones(data["N_missing"]),
            "alpha_sq": 500**2 * np.ones(data["N"]),
            "variance_sys": scale**2 * np.ones(data["N"]),
            "L_corr": np.eye(data["N"])
        }

        raise a
        op_params = model.optimize(data, init=init, iter=100000)

        raise a
        fit = model.sample(data, init=op_params, **sample_kwds)

        model.write(model_path, overwrite=True)

        # Keep the model.
        models[wg][parameter] = model


        # Plot some stuff.
        fig = model.plot_node_correlations()
        
        raise a
        fig.savefig("wg{}-{}-node-correlations.pdf".format(wg, parameter))

        fig = model.plot_node_uncertainty_with_snr(
            show_data_points=False, legend_kwds=dict(ncol=2))
        fig.savefig("wg{}-{}-node-uncertainties.pdf".format(wg, parameter))

        continue



        # Homogenise this parameter for all stars.
        model.homogenise_all_stars(update_database=True)

        continue        


# TODO:
# - VROT
# - VEL
# - FLAGS
# - XI

# Some notes:
# - If a node provides `mh` or `feh`, we can treat them as the same (see `scripts/setup_db.py`)
