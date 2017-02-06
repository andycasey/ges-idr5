
"""
Homogenise the results from "WG20" (WG10 HR15N)
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

CHECK_FOR_NODE_SETUPS = False

# Create a database object.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)

# Give new node names to all WG10 nodes on a per-setup basis?
if CHECK_FOR_NODE_SETUPS:

    distinct_nodes = database.retrieve_table(
        """select distinct on (name, node_id, setup) name, node_id, setup from results, nodes where results.node_id = nodes.id and nodes.wg = 20 and results.teff <> 'NaN'""")

    if not np.all(["-" in name for name in distinct_nodes["name"]]):
        logger.info("Creating new nodes!")

        for distinct_node in distinct_nodes:
            new_node_name = "{}-{}".format(distinct_node["name"].strip(), distinct_node["setup"].strip())
            new_node_id = database.create_or_retrieve_node_id(10, new_node_name)

            # Update all the results.
            logger.info("Updating to {}/{}".format(new_node_id, new_node_name))
            N = database.update(
                "UPDATE results SET node_id = '{}' WHERE node_id = '{}';".format(
                    new_node_id, distinct_node["node_id"]))
            
            print(N)

        database.connection.commit()

# Load the "benchmarks"
# Only use "benchmarks" with TEFF < 8000 K
benchmarks = Table.read("fits-templates/benchmarks/GES_iDR5_FGKM_Benchmarks_ARC_29092016.fits")
benchmarks = benchmarks[benchmarks["TEFF"] < 8000]
benchmarks["E_FEH"] = 0.10

model_paths = "homogenisation-wg{wg}-{parameter}.model"


wgs = (20, )

parameter_scales = OrderedDict([
    ("teff", 250),
    ("logg", 0.25),
    ("feh", 0.25),

])

sample_kwds = dict(chains=4, iter=5000)

finite = np.isfinite(benchmarks["TEFF"] * benchmarks["LOGG"] * benchmarks["FEH"])
benchmarks = benchmarks[finite]


models = {}
for wg in wgs:
    models[wg] = {}
    for parameter, scale in parameter_scales.items():

        model_path = model_paths.format(wg=wg, parameter=parameter)
        model_path = "homogenisation-wg{wg}-{parameter}-dual-poly-sys.model".format(
            wg=wg, parameter=parameter)

        if os.path.exists(model_path): 
            model = EnsembleModel.read(model_path, database,
                model_path="code/model/ensemble-model-dual-2.stan")

        else:
            model = EnsembleModel(database, wg, parameter, benchmarks,
                model_path="code/model/ensemble-model-dual-2.stan")
            data, metadata = model._prepare_data(
                default_sigma_calibrator=scale, minimum_node_estimates=1)

            init = {
                "truths": data["mu_calibrator"],
                "biases": np.zeros(data["N"]),
                "missing_estimates": np.random.uniform(
                    data["lower_bound"], data["upper_bound"], size=data["TM"]),
                "alpha_sq": np.mean([data["lower_alpha_sq"], data["upper_alpha_sq"]]) * np.ones(data["N"]),
                "vs_c": scale**2 * np.ones(data["N"]),
                #"vs_a": 1e-2 + np.zeros((data["N"], data["S"])).T,
                #"vs_b": 1e-2 + np.ones((data["N"], data["S"])).T,
                "vs_theta": 1e-2 + np.zeros((9, data["N"])),
                "L_corr": np.eye(data["N"])
            }

            init.update({"vs_tb{}".format(i): 1 for i in range(1, 9)})
            init.update({"vs_lb{}".format(i): 1 for i in range(1, 9)})
            init.update({"vs_fb{}".format(i): 1 for i in range(1, 9)})

            init.update({"vs_ta{}".format(i): 2 for i in range(1, 9)})
            init.update({"vs_la{}".format(i): 2 for i in range(1, 9)})
            init.update({"vs_fa{}".format(i): 2 for i in range(1, 9)})
            
            init["vs_tc7"] = 1 * np.ones(data["N"])
            init["vs_tc8"] = 1 * np.ones(data["N"])
            init["vs_tc9"] = 1 * np.ones(data["N"])


            init.update({"sv_tb{}".format(i): 1 for i in range(1, 9)})
            init.update({"sv_lb{}".format(i): 1 for i in range(1, 9)})
            init.update({"sv_fb{}".format(i): 1 for i in range(1, 9)})

            init.update({"sv_ta{}".format(i): 2 for i in range(1, 9)})
            init.update({"sv_la{}".format(i): 2 for i in range(1, 9)})
            init.update({"sv_fa{}".format(i): 2 for i in range(1, 9)})
            
            init["sv_tc7"] = 1 * np.ones(data["N"])
            init["sv_tc8"] = 1 * np.ones(data["N"])
            init["sv_tc9"] = 1 * np.ones(data["N"])

            print("Number of model parameters for {}: {}".format(parameter,
                sum([np.array(v).size for v in init.values()])))

            op_params = model.optimize(data, init=init, iter=100000)
            op_params["sv_ta2"] += 1e-3

            fit = model.sample(data, init=op_params, **sample_kwds)
            #print("Printing fit")
            #print(fit)
            fig = fit.plot(pars=("biases", ))
            fig.savefig("homogenisation-wg{}-{}.png".format(wg, parameter))

            model.write(model_path, 
                overwrite=True, __ignore_model_pars=("Sigma", "full_rank_estimates"))

        #models[wg][parameter] = model

        #model.homogenise_benchmark_stars(update_database=True)
        #model.homogenise_all_stars(update_database=True, autocommit=True)



# TODO:
# - VROT
# - VEL
# - FLAGS
# - XI

# Some notes:
# - If a node provides `mh` or `feh`, we can treat them as the same (see `scripts/setup_db.py`)
