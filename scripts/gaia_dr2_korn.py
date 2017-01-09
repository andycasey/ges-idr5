
""" Extract calibration/validation targets for Gaia. """

__author__ = "Andrew R. Casey (arc@ast.cam.ac.uk)"

import yaml
from code import GESDatabase

# Create a database object.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)

table = database.retrieve_table(
    """ SELECT  DISTINCT ON (s.cname) s.ra, s.dec, s.cname
          FROM  spectra AS s,
                wg_recommended_results AS r 
         WHERE  r.cname = s.cname 
           AND  s.is_blind_test = false
           AND  r.cname <> 'ssssssss-sssssss'
           AND  (
                    (s.setup IN ('U580', 'U520') AND s.snr > 50)
                OR  (s.setup IN ('HR10', 'HR21') AND s.snr > 50)
                )
           AND  r.wg IN (10, 11)
           AND  s.vel <> 'NaN'
           AND  s.e_vel <> 'NaN'
           AND  s.vel > -500 AND s.vel < 500
           AND  s.ges_type NOT LIKE 'AR%'
           AND  r.teff > 4000
           AND  r.teff < 6500
           AND  r.e_teff < 250 
           AND  r.nn_nodes_teff >= 2""")

print(len(table))

table.write("ges-dr5-validation-for-gaia-30112016.csv")
