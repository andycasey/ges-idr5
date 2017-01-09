
"""
What is up with sub-giant stars in NGC2420?
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

t2 = database.retrieve_table("SELECT DISTINCT ON (r.id) r.id, r.cname, r.wg, s.ra, s.dec, s.vel, s.e_vel, nn_nodes_teff, nn_nodes_logg, teff, logg, feh, e_teff, e_logg, e_feh from wg_recommended_results as r, spectra as s where r.cname = s.cname and s.ges_fld like 'NGC2420%' and r.wg = 11;")
t2.write("sandbox/ngc2420-homogenised.fits", overwrite=1)


t = database.retrieve_table("SELECT distinct on (r.id) r.id, r.cname, s.ra, s.dec, s.vel, s.e_vel, n.name, teff, logg, feh, passed_quality_control from results as r, nodes as n, spectra as s where r.node_id = n.id and r.cname = s.cname and s.ges_fld like 'NGC2420%' and n.wg = 11 and passed_quality_control;")
t.write("sandbox/ngc2420-nodes.fits")



raise a
