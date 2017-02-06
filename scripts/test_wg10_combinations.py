
"""
Test different mappings between WG10 nodes.
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


# Only use "benchmarks" with TEFF < 8000 K
benchmarks = Table.read("fits-templates/benchmarks/GES_iDR5_FGKM_Benchmarks_ARC_29092016.fits")
benchmarks = benchmarks[benchmarks["TEFF"] < 8000]
benchmarks["E_FEH"] = 0.10

vector_terms = {
    "teff":[
        (("teff", ), 1),
        (("logg", ), 1),
        (("feh", ), 1),
        (("teff", "logg", ), 1),
        (("teff", "feh", ), 1),
        (("logg", "feh", ), 1),
        (("teff", ), 2),
        (("logg", ), 2),
        (("feh", ), 2),
    ],
    "logg": [
        (("teff", ), 1),
        (("logg", ), 1),
        (("feh", ), 1),
        (("teff", "logg", ), 1),
        (("teff", "feh", ), 1),
        (("logg", "feh", ), 1),
        (("teff", ), 2),
        (("logg", ), 2),
        (("feh", ), 2)
    ],
    "feh": [
        (("feh", ), 1),
    ]

}


models = {}
for parameter in ("teff", "logg", "feh"):


    model = CombinedModel(database, 10, parameter, benchmarks)

    data, metadata = model._prepare_mapped_data(
        vector_terms=vector_terms[parameter],
        reference_node_id=database.retrieve_node_id(10, "Lumba-HR10|HR21"),
        minimum_node_estimates=3,
        sql_constraint_for_mapping_query=\
            "r.snr > 10 and r.teff > 4000 and r.teff < 6250 and "
            "r.logg > 0 and r.logg < 5 and r.feh > -1.5 and r.feh < +0.25"
            " and passed_quality_control",
        sql_constraint_for_data_query=None)

    models[parameter] = model


# Get the mapped data for each node.
for node_name in metadata["node_names"]:

    if "Nice" not in node_name: continue

    node_id = database.retrieve_node_id(10, node_name)

    # Apply the coefficients.
    mapped_values = {}
    for parameter, model in models.items():
        mapped_data = model._apply_mapping_coefficients(sql_query="""
            SELECT cname, node_id, teff, e_teff, logg, e_logg, feh, e_feh, snr
            FROM results
            WHERE node_id = '{}' and passed_quality_control;
            """.format(node_id))

        mapped_values[parameter] = mapped_data[parameter]


    data = database.retrieve_table("""
        SELECT cname, node_id, teff, e_teff, logg, e_logg, feh, e_feh, snr
        FROM results
        WHERE node_id = '{}' and passed_quality_control;
        """.format(node_id))

    for parameter, values in mapped_values.items():
        data[parameter] = values

    data.write(
        "node-results/ges-dr5-wg10-MAPPED-{}.fits".format(node_name.split("-")[0]),
        overwrite=True)

    print("Done with {}".format(node_name))
