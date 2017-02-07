
"""
A script to homogenise Gaia-ESO Survey data using a combined model,
one that maps node results onto a common scale before inferring their
relative contributions and correlations using a noise model.
"""

import yaml
import logging
import numpy as np
import os
from astropy.table import Table
from collections import OrderedDict


from code import GESDatabase
from code.model.combined import CombinedModel


# Initialize logging.
logger = logging.getLogger("ges")

# Create a database object.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)

# Enable autocommit mode so that we can do homogenisation in parallel
database.connection.autocommit = True

# Load the "benchmarks"
# Only use "benchmarks" with TEFF < 8000 K
benchmarks = Table.read("fits-templates/benchmarks/GES_iDR5_FGKM_Benchmarks_ARC_29092016.fits")
benchmarks = benchmarks[benchmarks["TEFF"] < 8000]
benchmarks["E_FEH"] = 0.10

model_paths = "homogenisation-wg{wg}-{parameter}-{setup}-combined-wrt-MaxPlanck.model"

wgs = (10, )
parameter_scales = OrderedDict([
#    ("teff", 250),
#    ("logg", 0.25),
    ("feh", 0.25)
])


vector_terms = {
    "teff":[
        (("teff", ), 1),
        (("logg", ), 1),
        (("feh", ), 1),
    ],
    "logg": [
        (("teff", ), 1),
        (("logg", ), 1),
        (("feh", ), 1),
    ],
    "feh": [
        (("feh", ), 1),
    ]

}

reference_node_id = database.retrieve_node_id(10, "MaxPlanck-HR21")


sample_kwds = dict(chains=2, iter=4000)

finite = np.isfinite(benchmarks["TEFF"] * benchmarks["LOGG"] * benchmarks["FEH"])
benchmarks = benchmarks[finite]

for wg in wgs:
    for setup in ("HR10", "HR21", "HR10|HR21"):

        for parameter, default_sigma_calibrator in parameter_scales.items():

            model_path = model_paths.format(
                wg=wg, parameter=parameter, setup=setup.replace("|", "_"))
            
            if os.path.exists(model_path): 
                model = CombinedModel.read(model_path, database)
                
            else:
                model = CombinedModel(database, 10, parameter, benchmarks)
                data, metadata = model._prepare_mapped_data(
                    default_sigma_calibrator=default_sigma_calibrator, 
                    minimum_node_estimates=1,
                    vector_terms=vector_terms[parameter],
                    reference_node_id=reference_node_id,
                    sql_constraint_for_mapping_query="(r.teff < 6000 or name not like 'Nice%') and r.snr > 50 and r.teff > 4000 and r.teff < 6500 and "
                        "r.logg > 0 and r.logg < 5 and r.feh > -1 and r.feh < +0.25"
                        " and passed_quality_control and (trim(name) like '%-{setup}' or node_id = {reference_node_id})".format(
                            setup=setup, reference_node_id=reference_node_id),
                    sql_constraint_for_data_query="TRIM(name) LIKE '%-{}'".format(setup),
                    match_by="cname")


                init = {
                    "truths": data["mu_calibrator"],
                    "biases": np.zeros(data["N"]),
                    "missing_estimates": np.random.uniform(
                        data["lower_bound"], data["upper_bound"], size=data["TM"]),
                    "alpha_sq": np.mean([data["lower_alpha_sq"], data["upper_alpha_sq"]]) * np.ones(data["N"]),
                    "systematic_variance": default_sigma_calibrator**2 * np.ones(data["N"]),
                    "L_corr": np.eye(data["N"])
                }

                op_params = model.optimize(data, init=init, iter=100000)

                fit = model.sample(data, init=op_params, **sample_kwds)

                model.write(model_path, 
                    overwrite=True, __ignore_model_pars=("Sigma", "full_rank_estimates"))

            model.homogenise_stars_matching_query("""
                SELECT DISTINCT ON (cname) cname
                FROM  results, nodes
                WHERE results.node_id = nodes.id
                  AND nodes.wg = '{wg}'
                  AND TRIM(nodes.name) LIKE '%-{setup}'
                  AND {parameter} <> 'NaN'
                  AND passed_quality_control
                """.format(wg=wg, setup=setup, parameter=parameter),
                sql_constraint="TRIM(nodes.name) LIKE '%-{}'".format(setup),
                metadata={"setup": setup})


        # ANDY: DO HR10 TEFF WG10 -- it is not complete.