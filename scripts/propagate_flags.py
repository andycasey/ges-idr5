#!/usr/bin/python

"""
Propagate relevant flag information from one node to others.
"""

import logging
import yaml

from code import GESDatabase

# Connect to database.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)

# Create a database object.
database = GESDatabase(**credentials)

logger = logging.getLogger("ges")


with open("flags.yaml", "r") as fp:
    qc_flags = yaml.load(fp)


# Clear any previous propagations before starting.
'''
logger.info("Clearing previous propagations and setting all to have passed_quality_control = True")
database.update(
    """ UPDATE  results
           SET  propagated_tech_from_result_id = null,
                propagated_peculi_from_result_id = null,
                propagated_remark_from_result_id = null,
                propagated_tech = '',
                propagated_peculi = '',
                propagated_remark = '',
                passed_quality_control = true;""")
'''

N_propagations = {}
N_marked_as_poor_quality = {}

# PROPAGATE BY SPECTRUM

for flag in qc_flags["propagate_flags_by_spectrum"]:

    N_propagations.setdefault(flag, 0)

    affected = database.retrieve_table(
        """ SELECT  id, TRIM(filename) AS filename, TRIM(tech) as tech
            FROM    results
            WHERE   tech LIKE '%{}-%'
        """.format(flag))
    
    if affected is None:
        continue

    for row in affected:
    
        # Each row can have multiple TECH flags, so first identify the TECH flag
        # that PostgreSQL matched on.
        tech_flags = row["tech"].strip().split("|")
        for tech_flag in tech_flags:
            if "{}-".format(flag) in tech_flag:
                matched_tech_flag = tech_flag
                break

        else:
            raise ValueError(
                "cannot identify tech flag {} from the SQL match: {}".format(
                    flag, row["tech"].strip()))


        # Update other results using the same filename(s).
        filenames = row["filename"].strip().split("|")
        for j, filename in enumerate(filenames):
            if not filename: continue

            N = database.update(
                """ UPDATE  results
                    SET     propagated_tech_from_result_id = '{}',
                            propagated_tech = '{}',
                            passed_quality_control = false
                    WHERE   filename LIKE '%{}%'
                      AND   passed_quality_control = true
                """.format(
                    int(row["id"]), matched_tech_flag, filename))
            
            N_propagations[flag] += N

            if N > 0:
                logger.info("Propagated ({}/{}/{}) to {} other entries".format(
                    row["id"], matched_tech_flag, filename, N))    


        database.connection.commit()


# PROPAGATE BY CNAME

for flag in qc_flags["propagate_flags_by_cname"]:

    N_propagations.setdefault(flag, 0)

    affected = database.retrieve_table(
        """ SELECT  id, cname, TRIM(tech) as tech
              FROM  results
             WHERE  tech LIKE '%{}-%'
        """.format(flag))

    if affected is None:
        continue

    for row in affected:
        # Each row can have multiple TECH flags, so first identify the TECH flag
        # that PostgreSQL matched on.
        tech_flags = row["tech"].strip().split("|")
        for tech_flag in tech_flags:
            if "{}-".format(flag) in tech_flag:
                matched_tech_flag = tech_flag
                break

        else:
            raise ValueError(
                "cannot identify tech flag {} from the SQL match: {}".format(
                    flag, row["tech"].strip()))

        # Update other results matching this CNAME.
        N = database.update(
            """ UPDATE  results
                   SET  propagated_tech_from_result_id = '{}',
                        propagated_tech = '{}',
                        passed_quality_control = false
                 WHERE  cname = '{}'
                   AND  passed_quality_control = true;
            """.format(int(row["id"]), matched_tech_flag, row["cname"]))
        N_propagations[flag] += N

        if N > 0:
            logger.info("Propagated ({}/{}/{}) to {} other entries".format(
                row["id"], matched_tech_flag, row["cname"], N))

        database.connection.commit()


# NODE-SPECIFIC FLAGS

for node_specific_flag in qc_flags["node_specific_flags"]:

    N = database.update(
        """ UPDATE  results
               SET  passed_quality_control = false
             WHERE  tech LIKE '%{}-%'
               AND  passed_quality_control = true;
        """.format(node_specific_flag))
    N_marked_as_poor_quality[node_specific_flag] = N

    if N > 0:
        logger.info("Marked {} results as poor quality due to matching flag {}"\
            .format(N, node_specific_flag))

database.connection.commit()