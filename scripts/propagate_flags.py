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
logger.info("Clearing previous propagations and setting all to have passed_quality_control = True")
database.update(
    """ UPDATE  results
           SET  propagated_tech_from_result_id = null,
                propagated_peculi_from_result_id = null,
                propagated_remark_from_result_id = null,
                propagated_tech = '',
                propagated_peculi = '',
                propagated_remark = '',
                passed_quality_control = true
         WHERE  passed_quality_control = false;""")

database.connection.commit()

# Identify spurious spectra and mark them as such.
N_peculiar_spectra = {}
peculiar_spectra_kwds = dict(
    and_or="and", sigma_discrepant=3, teff_discrepant=250, logg_discrepant=0.25)

for wg in (10, 11, 12, 13):

    logger.info("Querying for peculiar spectra in WG{}".format(wg))

    kwds = dict(wg=wg)
    kwds.update(peculiar_spectra_kwds)

    peculiar_spectrum_query = """
        with t4 as (
        select id, cname, filename, avg_filename_teff, avg_cname_teff, stddev_cname_teff, abs((avg_cname_teff - avg_filename_teff)/(0.00001 + stddev_cname_teff)) as abs_sigma_teff_discrepant, avg_filename_logg, avg_cname_logg, stddev_cname_logg, abs((avg_cname_logg - avg_filename_logg)/(0.00001 + stddev_cname_logg)) as abs_sigma_logg_discrepant FROM (with ar as (select distinct on (filename) id, cname, trim(filename) as filename, avg(teff) over w as avg_filename_teff, avg(logg) over w as avg_filename_logg from (with n as (select id from nodes where wg = {wg}) select distinct on (r.filename, r.node_id) r.id, r.cname, trim(r.filename) as filename, r.node_id, teff, logg from n, results as r where r.node_id = n.id and teff <> 'NaN' or logg <> 'NaN') t window w as (partition by filename)) select ar.id, ar.cname, ar.filename, ar.avg_filename_teff, avg(avg_filename_teff) over w2 as avg_cname_teff, stddev(avg_filename_teff) over w2 as stddev_cname_teff, ar.avg_filename_logg, avg(avg_filename_logg) over w2 as avg_cname_logg, stddev(avg_filename_logg) over w2 as stddev_cname_logg FROM ar window w2 as (partition by cname)) t3)
        select * from t4 where (t4.abs_sigma_teff_discrepant > {sigma_discrepant} and t4.abs_sigma_teff_discrepant <> 'NaN' and abs(t4.avg_cname_teff - avg_filename_teff) >= {teff_discrepant}) {and_or} (t4.abs_sigma_logg_discrepant > {sigma_discrepant} and abs(t4.avg_cname_logg - t4.avg_filename_logg) >= {logg_discrepant} and t4.abs_sigma_logg_discrepant <> 'NaN') order by cname asc;""".format(
            **kwds)

    peculiar_spectra = database.retrieve_table(peculiar_spectrum_query)
    N_peculiar_spectra[wg] = len(peculiar_spectra)

    for row in peculiar_spectra:

        filenames = row["filename"].strip().split("|")
        logger.info("Propagating {}/{}/{}".format(
            row["id"], row["cname"], row["filename"]))

        n = 0
        for filename in filenames:
            n += database.update(
                """ UPDATE results
                       SET propagated_tech = '10106-{}-00-00-A',
                           propagated_tech_from_result_id = '{}',
                           passed_quality_control = false
                     WHERE filename LIKE '%{}%';""".format(
                        wg, int(row["id"]), filename))

        if n > 0:
            logger.info("--> affected {} results".format(n))

database.connection.commit()


# PROPAGATE BY SPECTRUM
def propagate_by_spectrum(flag, constraint=None):
    
    constraint_str = "" if constraint is None else " AND {}".format(constraint)
    affected = database.retrieve_table(
        """ SELECT  id, TRIM(filename) AS filename, TRIM(tech) as tech
            FROM    results
            WHERE   tech LIKE '%{}-%' {}
        """.format(flag, constraint_str))
    
    if affected is None:
        return 0

    N = 0
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

            n = database.update(
                """ UPDATE  results
                    SET     propagated_tech_from_result_id = '{}',
                            propagated_tech = '{}',
                            passed_quality_control = false
                    WHERE   filename LIKE '%{}%'
                      AND   passed_quality_control = true
                """.format(
                    int(row["id"]), matched_tech_flag, filename))
            
            N += n

            if n > 0:
                logger.info("Propagated ({}/{}/{}) to {} other entries".format(
                    row["id"], matched_tech_flag, filename, n))    

        database.connection.commit()

    return N


# PROPAGATE BY SPECTRUM
N_propagations = {}

for key, value in qc_flags["propagate_flags_by_spectrum"].items():
    if key == "no_constraint":
        for flag in value:
            N_propagations.setdefault(flag, 0)
            N_propagations[flag] += propagate_by_spectrum(flag, None)
    
    else:
        flag = key
        constraint = value.get("constraint", None)
        N_propagations.setdefault(flag, 0)
        N_propagations[flag] += propagate_by_spectrum(flag, constraint)



# PROPAGATE BY CNAME
def propagate_by_cname(flag, constraint=None):

    constraint_str = "" if constraint is None else " AND {}".format(constraint)
    affected = database.retrieve_table(
        """ SELECT  id, cname, TRIM(tech) as tech
              FROM  results
             WHERE  tech LIKE '%{}-%' {}
        """.format(flag, constraint_str))

    if affected is None:
        return 0

    N = 0
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
        n = database.update(
            """ UPDATE  results
                   SET  propagated_tech_from_result_id = '{}',
                        propagated_tech = '{}',
                        passed_quality_control = false
                 WHERE  cname = '{}'
                   AND  passed_quality_control = true;
            """.format(int(row["id"]), matched_tech_flag, row["cname"]))
        N += n

        if n > 0:
            logger.info("Propagated ({}/{}/{}) to {} other entries".format(
                row["id"], matched_tech_flag, row["cname"], n))

        database.connection.commit()

    return N


# PROPAGATE BY CNAME
for key, value in qc_flags["propagate_flags_by_cname"].items():

    if key == "no_constraint":
        for flag in value:
            N_propagations.setdefault(flag, 0)
            N_propagations[flag] += propagate_by_cname(flag, None)

    else:
        flag = key
        constraint = value.get("constraint", None)
        N_propagations.setdefault(flag, 0)
        N_propagations[flag] += propagate_by_cname(flag, constraint)


# NODE-SPECIFIC FLAGS
N_marked_as_poor_quality = {}
for key, value in qc_flags["node_specific_flags"].items():

    if key == "no_constraint":
        for flag in value:
            N_marked_as_poor_quality.setdefault(flag, 0)
            N = database.update(
                """ UPDATE  results
                       SET  passed_quality_control = false
                     WHERE  tech LIKE '%{}-%'
                       AND  passed_quality_control = true;
                """.format(flag))
            N_marked_as_poor_quality[flag] += N

            if N > 0:
                logger.info(
                    "Marked {} results as poor quality due to matching flag {}"\
                    .format(N, flag))
    else:
        flag = key
        constraint = value.get("constraint", None)
        constraint_str = "" if constraint is None else " AND {}".format(constraint)
        N_marked_as_poor_quality.setdefault(flag, 0)

        N = database.update(
            """ UPDATE  results
                   SET  passed_quality_control = false
                 WHERE  tech LIKE '%{}-%'
                   AND  passed_quality_control = true
                   {};
            """.format(flag, constraint_str))
        N_marked_as_poor_quality[flag] += N

        if N > 0:
            logger.info(
                "Marked {} results as poor quality due to matching flag {} "\
                "and constraint {}".format(N, flag, constraint))

    database.connection.commit()