
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


# Assign the nodes new names first.
old_nodes = database.retrieve_table("""
    SELECT DISTINCT ON (nodes.name, setup) nodes.id, nodes.name, setup
    FROM results, nodes
    WHERE results.node_id = nodes.id
      AND nodes.wg = 10
      AND results.teff <> 'NaN'""")

for old_node in old_nodes:

    if "-" in old_node["name"]:
        print("Skipping {}".format(old_node["name"]))
        continue

    new_node_name = "{name}-{setup}".format(
        name=old_node["name"], setup=old_node["setup"].strip())

    new_node_id = database.create_or_retrieve_node_id(10, new_node_name)

    N = database.update("""
        UPDATE results 
        SET node_id = %s
        WHERE node_id = %s
          AND TRIM(setup) = %s;
        """, (new_node_id, old_node["id"], old_node["setup"].strip()))

    print(old_node["id"], old_node["name"], old_node["setup"].strip(), N)


# Only use "benchmarks" with TEFF < 8000 K
benchmarks = Table.read("fits-templates/benchmarks/GES_iDR5_FGKM_Benchmarks_ARC_29092016.fits")
benchmarks = benchmarks[benchmarks["TEFF"] < 8000]
benchmarks["E_FEH"] = 0.10

# Do everything except Nice first.
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

reference_nodes = {
    "HR10": ("Lumba-HR10|HR21", "cname"),
    "HR21": ("Lumba-HR10|HR21", "cname"),
    "HR10|HR21": ("Lumba-HR10|HR21", "filename"),
}

for setup, (reference_node, match_by) in reference_nodes.items():

    reference_node_id = database.retrieve_node_id(10, reference_node)

    models = {}
    for parameter in ("teff", "logg", "feh"):

        model = CombinedModel(database, 10, parameter, benchmarks)
        coeffs = model._get_mapping_coefficients(
            vector_terms=vector_terms[parameter], 
            reference_node_id=reference_node_id,
            sql_constraint="r.snr > 10 and r.teff > 4000 and r.teff < 6250 and "
                "r.logg > 0 and r.logg < 5 and r.feh > -1.5 and r.feh < +0.25"
                " and passed_quality_control and (trim(name) like '%-{}' or node_id = {})".format(setup, reference_node_id),
            match_by=match_by)

        models[parameter] = (model, coeffs)

    # Get the mapped data for each node.
    node_ids = models.values()[0][-1].keys()

    for node_id in node_ids:

        node_name = database.retrieve(
            "SELECT name FROM nodes WHERE id = {}".format(node_id))[0][0]

        if node_name.startswith("MaxPlanck"):
            for parameter in models.keys():
                # Don't do any mappings for MaxPlanck
                models[parameter][1][node_id] = models[parameter][1][reference_node_id]

        # Apply the coefficients.
        mapped_values = {}
        for parameter, (model, coeffs) in models.items():
            mapped_data = model._apply_mapping_coefficients(sql_query="""
                SELECT cname, node_id, teff, e_teff, logg, e_logg, feh, e_feh, snr
                FROM results
                WHERE node_id = '{}' and passed_quality_control;
                """.format(node_id),
                coefficients=coeffs, vector_terms=vector_terms[parameter])

            mapped_values[parameter] = mapped_data[parameter]

        data = database.retrieve_table("""
            SELECT id, cname, node_id, teff, e_teff, logg, e_logg, feh, e_feh, snr
            FROM results
            WHERE node_id = '{}' and passed_quality_control;
            """.format(node_id))

        # Update the database with the mapped values.
        for j, result_id in enumerate(data["id"]):

            database.execute("""
                UPDATE results
                SET teff = %s,
                    logg = %s,
                    feh = %s,
                    passed_quality_control = %s
                WHERE id = %s
                """, (data["teff"][j], data["logg"][j], data["feh"][j],
                    False if (node_name.startswith("IAC") and data["logg"][j] > 3.5) else True,
                    result_id))
        
        for parameter, values in mapped_values.items():
            data[parameter] = values

        data.write(
            "node-results/ges-dr5-wg10-MAPPED-{}.fits".format(node_name.replace("|", "-")),
            overwrite=True)

        print("Done with {}".format(node_name))


all_wg_10_node_ids = list(database.retrieve_table(
    "SELECT id FROM nodes WHERE wg = 10")["id"])

database.execute("""
    UPDATE results
    SET passed_quality_control = false
    WHERE node_id in %s
    AND (
        (teff < 4000) or
        (teff > 6500) or
        (logg > 5) or
        (logg < 0) or
        (feh > 0.5) or
        (feh < -2.5))
    """, (tuple(all_wg_10_node_ids), ))

# Mark Nice HR21 results with teff < 4500 and logg > 4 as bad.
node_id = database.retrieve_node_id(10, "Nice-HR21")
database.update(
    """ UPDATE  results
        SET     passed_quality_control = false
        WHERE   node_id = %s
          AND   logg > 4
          AND   teff < 4500
    """, (node_id, ))

database.connection.commit()
