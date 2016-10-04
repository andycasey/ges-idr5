
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

wgs = (11, )
parameter_scales = OrderedDict([
    ("teff", 250),
    ("logg", 0.25),
    ("feh", 0.25)
])

sample_kwds = dict(chains=4, iter=2000)

benchmarks = benchmarks[(benchmarks["FEH"] > -2) * (benchmarks["TEFF"] > 4000)]

models = {}
for wg in wgs:
    models[wg] = {}
    for parameter, scale in parameter_scales.items():

        model_path = model_paths.format(wg=wg, parameter=parameter)

        if wg == 11 and parameter == "teff": continue

        '''
        model = SingleParameterEnsembleModel.read(model_path, database)
        model._prepare_data()


        # SHow heaps of shit.
        data = database.retrieve_table(
            """ SELECT DISTINCT ON (node_id, filename) r.id, node_id, nodes.name, spectra.ges_fld, teff, r.snr, TRIM(r.filename) as filename
                FROM results as r, nodes, spectra
                WHERE trim(spectra.ges_type) like 'GE_SD_B%'
                  AND spectra.cname = r.cname
                  AND r.node_id = nodes.id
                  AND nodes.wg = '11'
                  AND passed_quality_control
                  AND teff <> 'NaN'
            """)

        data = data.group_by("node_id")

        for group in data.groups:
            fig, ax = plt.subplots()

            diffs = []
            keep = np.zeros(len(group), dtype=bool)
            ges_flds = np.array([each.strip() for each in benchmarks["GES_FLD"]])
            for i, row in enumerate(group):
                match = ges_flds == row["ges_fld"].strip()
                if any(match):
                    diffs.append(row["teff"] - benchmarks["TEFF"][match])
                    keep[i] = True

            diffs = np.abs(diffs)



            ax.scatter(group["snr"][keep], diffs, label=group["name"][0], alpha=0.5)

            match = np.where(group["node_id"][0] == model._metadata["node_ids"])[0][0]

            diffs = np.array([model._data["mu_calibrator"][index - 1]  for index in model._data["calibrator_index"]]) \
                - model._data["estimates"][:, match]

            ax.scatter(model._data["snr_spectrum"], np.abs(diffs), facecolor="r", alpha=0.5)





            wtf = (diffs.flatten() > 1000) * (group["snr"][keep] > 100)

            print(data[keep][wtf])

            ax.set_title(group["name"][0])


        raise a
        '''


        if os.path.exists(model_path):
            model = SingleParameterEnsembleModel.read(model_path, database)

        else:
            model = SingleParameterEnsembleModel(database, wg, parameter, benchmarks)
            data, metadata = model._prepare_data(default_sigma_calibrator=scale)

            init = {
                "truths": data["mu_calibrator"],
                "var_sys_estimator": (scale/5.)**2 * np.ones(data["N_estimators"]),
                "alpha": 750. * np.ones(data["N_estimators"]),
                "rho_estimators": np.zeros(metadata["N_pairwise_estimators"]),
                "c0_estimators": np.zeros(data["N_estimators"])
            }

            op_params = model.optimize(data, init=init)

            fit = model.sample(data, init=op_params, **sample_kwds)

            model.write(model_path, overwrite=True)

        fig = model.node_uncertainty_with_snr(show_data_points=False)
        fig.savefig("wg{}-{}-node-uncertainty.pdf".format(wg, parameter))
        continue

        # Homogenise this parameter for all stars.
        model.homogenise_all_stars(update_database=True)

        # Keep the model.
        models[wg][parameter] = model

# TODO:
# - VROT
# - VEL
# - FLAGS
# - XI

# Some notes:
# - If a node provides `mh` or `feh`, we can treat them as the same (see `scripts/setup_db.py`)
