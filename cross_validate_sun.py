
"""
Cross-validate ensemble model on the sun.
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


benchmarks["E_FEH"] = 0.2

bms_sun_excluded = benchmarks[benchmarks["TEFF"] < 8000]
keep = np.array([e.strip() != "SUN" for e in bms_sun_excluded["GES_FLD"]])
assert len(keep) - sum(keep) > 0
bms_sun_excluded = bms_sun_excluded[keep]


cv_sun_model = NewEnsembleModel(database, 11, bms_sun_excluded)

parameters = {
    "logg": 0.25**2,
    "teff": 100**2,
    "feh": 0.20**2
}

saved_fits = {}
posteriors = {}

for parameter, var_intrinsic in parameters.items():

    if parameter in posteriors: continue

    data, metadata = cv_sun_model._prepare_data(parameter)

    init = {
        "truths": data["mu_calibrator"],
        "var_intrinsic": var_intrinsic**2,
        "var_sys_estimator": var_intrinsic**2 * np.ones(data["N_estimators"]),
        "alpha_sq": 1000 * np.ones(data["N_estimators"]),
        "rho_estimators": np.zeros(metadata["N_pairwise_estimators"]),
        "c0_estimators": np.zeros(data["N_estimators"])
    }

    op_params = cv_sun_model.optimize(data, init=init)

    samples = cv_sun_model.sample(data, init=op_params, chains=6, iter=2000)

    saved_fits[parameter] = samples

    posteriors[parameter] = homogenise_survey_measurements(
        "ssssssss-sssssss", 11, parameter, samples, database)


response = ""
for k, v in posteriors.items():
    
    p = np.percentile(v, [16, 50, 84])
    response += "{0}: {1:.2f} ({2:.2f}, {3:.2f})\n".format(
        k, p[1], p[0] - p[1], p[2] - p[1])

print(response)
with open("cv_solar.log", "w") as fp:
    fp.write(response)


raise a




# 18 Sco:
homogenise_survey_measurements(
    "16153746-0822162", 11, "logg", samples, database)







homogenise_survey_measurements(
    "14153967+1910567", 11, "logg", samples, database)