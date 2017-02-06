
"""
What is up with MS stars with Teff > 6000 in NGC 6633??
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

# Create a database object.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)

t2 = database.retrieve_table("SELECT DISTINCT ON (r.id) r.id, r.cname, r.wg, s.ra, s.dec, s.vel, s.e_vel, teff, logg, feh, e_teff, e_logg, e_feh from wg_recommended_results as r, spectra as s where r.cname = s.cname and s.ges_fld like 'NGC6633%' and r.wg = 11;")
t2.write("sandbox/ngc6633-homogenised.fits", overwrite=1)


t = database.retrieve_table("SELECT distinct on (r.id) r.id, r.cname, s.ra, s.dec, s.vel, s.e_vel, n.name, teff, logg, feh, passed_quality_control from results as r, nodes as n, spectra as s where r.node_id = n.id and r.cname = s.cname and s.ges_fld like 'NGC6633%' and n.wg = 11 and passed_quality_control;")
t.write("sandbox/ngc6633-nodes.fits")



# Now re-run without Vilnius.
vilnius_node_id = database.retrieve_node_id(11, "Vilnius")
database.execute(
    """ UPDATE results SET passed_quality_control = false
        WHERE id IN (
            SELECT DISTINCT ON (r.id) r.id
            FROM results AS r, spectra AS s
            WHERE r.cname = s.cname
              AND r.node_id = '{0:.0f}'
              AND s.ges_fld LIKE 'NGC6633%'
              AND r.passed_quality_control
            );""".format(vilnius_node_id))

raise a
