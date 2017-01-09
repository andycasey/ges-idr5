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

raise a


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
def propagate_by_spectrum(flag, constraint=None, commit=False):
    
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

    if commit:
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

database.connection.commit()


# PROPAGATE BY CNAME
def propagate_by_cname(flag, constraint=None, commit=False):

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

    if commit:
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

database.connection.commit()


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


# Flag *CLEARLY* erroneous problems from the WG11/IACAIP node. Some of these will
node_id = database.retrieve_node_id(11, "IACAIP")

N = database.update(
    """ UPDATE  results
           SET  passed_quality_control = false,
                propagated_tech = '10308-11-00-00-A'
         WHERE  results.node_id = '{}'
           AND  passed_quality_control
           AND  (
                    (teff <> 'NaN' AND (teff <= 3000 OR teff >= 8000))
                    OR
                    (logg <> 'NaN' AND (logg <= 0 OR logg >= 5))
                    OR
                    (feh <> 'NaN' AND (feh <= -3 OR feh >= 1))
                );
    """.format(node_id))
logger.info("Marked {} WG11/IACAIP results as out-of-grid.")

raise Verify("Check that no nodes have 10105 propagated unles sthey are node = 8")
raise Verify("Remove bad OACT results")
"""
Code used:

update results set propagated_tech_from_result_id = null, propagated_tech = '', passed_quality_control = true where propagated_tech = '10105-11-08-00-A' and node_id != 8 and node_id != 10;

3084 updated.
"""

# OACT is damaging.

# node_id = 11 is WG11/OACT
""" UPDATE  results
    SET     passed_quality_control = false
    WHERE   node_id = 11
      AND   feh <= -0.5;
"""
# 672 rows affected


# Since we can either just use a full record (e.g., TEFF/LOGG/FEH) or none of it,
# we have to do the uncomfortable thing and just set all OACT logg results to NaN
""" UPDATE  results
    SET     logg = 'NaN'
    WHERE   node_id = 11;
"""
# 6090 rows affected


# Remove archival results for 5927 and M12
"""UPDATE results SET passed_quality_control = false
WHERE cname IN (
    SELECT distinct on (r.cname) r.cname FROM nodes as n, results as r, spectra as s where r.node_id = n.id and r.cname = s.cname and (s.ges_fld like 'M12%' or s.ges_fld like 'NGC5927%') and n.wg = 11 and ges_type like 'AR_%'
    );
"""
#
 TODO:


# Lumba cluster stars are a bit screwy
# December 1, 2016
#select distinct on (delta_teff, cname) r1.passed_quality_control as lumba_qc, r2.passed_quality_control as nice_qc, r1.teff - r2.teff as delta_teff, r1.logg - r2.logg as delta_logg, r1.feh - r2.feh as delta_feh, r1.cname, s.ges_type, s.ges_fld from results as r1, results as r2, spectra as s where r1.cname = r2.cname and r1.node_id = 2 and r2.node_id = 10 and r1.tech like '10210%' and r1.cname = s.cname and r1.teff <> 'NaN' and r2.teff <> 'NaN' and ges_type not like '%_SD_B%' and ges_type like 'GE_CL%' order by delta_teff, cname desc;

"""
UPDATE results set passed_quality_control = false 
WHERE node_id = 2 and cname in (select r1.cname from results as r1, results as r2, spectra as s
    where r1.cname = r2.cname and r1.node_id = 2 and r2.node_id = 10 and s.cname = r1.cname and s.ges_type not like '%_SD_B%' and ges_type like 'GE_CL%' and abs(r1.teff - r2.teff) > 1000 and r1.teff <> 'NaN' and r2.teff <> 'NaN');
"""

# Lumba bad convergence for M15 stars
"""
DELETE FROM wg_recommended_results where cname in ('21300431+1210561', '21300738+1210330');
UPDATE results SET passed_quality_control = false WHERE  cname in ('21300431+1210561', '21300738+1210330') and node_id = 2;
"""

# Remove OACT and Vilnius for Rup134 to see if that improves clump.
"""
ges_idr5=# select distinct on (r.id) r.id from results as r, spectra as s, nodes as n where r.cname = s.cname and r.node_id = n.id and (n.id = 6 or n.id = 11) and r.passed_quality_control and s.ges_fld like 'Rup134%';
   id   
--------
 510895
 510896
 510897
 510898
 510899
 510900
 510901
 510902
 510903
 510904
 510905
 510906
 510907
 510908
 510909
 510910
 510911
 510912
 510913
 510914
 510915
 510916
 510917
 510918
 510919
 510920
 510921
 510922
 510923
 510924
 510925
 510926
 510927
 510928
 510929
 510930
 510931
 510932
 516985
 516986
 516987
 516988
 516989
 516990
 516991
 516992
 516993
 516994
 516995
 516996
 516997
 516998
 516999
 517000
 517001
 517002
 517003
 517004
 517005
 517006
 517007
 517009
 517010
 517011
 517012
 517013
 517014
 517015
 517016
 517017
 517018
 517019
 517020
 517021
 517022
(75 rows)
"""
"""
update results set passed_quality_control = false where id in (select distinct on (r.id) r.id from results as r, spectra as s, nodes as n where r.cname = s.cname and r.node_id = n.id and (n.id = 6 or n.id = 11) and r.passed_quality_control and s.ges_fld like 'Rup134%');
"""


# Identify outliers in the cluster measurements.
#select distinct on (r.id, abs_delta_teff) r.id, r.cname, n.name, trim(r.tech) as tech, s.ges_fld, s.ges_type, r.teff, wgr.teff as recommended_teff, wgr.e_teff as recommended_e_teff, wgr.nn_nodes_teff, wgr.nn_spectra_teff, abs(wgr.teff - r.teff) as abs_delta_teff from results as r, nodes as n, spectra as s, wg_recommended_results as wgr where s.cname = r.cname and r.cname = wgr.cname and abs(wgr.teff - r.teff) > 500 and wgr.wg = 11 and n.id = r.node_id and n.wg = 11 and r.teff <> 'NaN' and wgr.teff <> 'NaN' and r.passed_quality_control and (s.ges_type like '%_OC%' or s.ges_type like '%_GC%' or s.ges_type like '%_CL%') order by abs_delta_teff desc, r.id asc;
# --> saved to cluster-outlier.fits



database.connection.commit()