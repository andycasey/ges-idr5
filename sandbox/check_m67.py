
"""
What is up with M67?
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

# Setting these as not passing quality control to test result:

"""
ges_idr5=# select distinct on (r.id) r.id, n.name, r.teff, r.logg, r.cname from results as r, nodes as n, spectra as s where r.cname = s.cname and r.node_id = n.id and n.wg = 11 and (n.name like 'Vilnius%' or n.name like 'Nice%') and s.ges_fld like 'M67%' and r.passed_quality_control and r.teff < 5250 and r.logg < 3.85;
   id   |  name   |       teff       |        logg        |      cname       
--------+---------+------------------+--------------------+------------------
 504050 | Nice    | 4876.60009765625 | 3.5920000076293945 | 08510838+1147121
 504060 | Nice    | 4850.39990234375 | 3.5160000324249268 | 08513045+1148582
 504062 | Nice    |  4878.2001953125 | 3.5399999618530273 | 08513577+1153347
 504066 | Nice    |           4745.5 | 3.1421999931335449 | 08514507+1147459
 516230 | Vilnius |           5005.0 | 3.7000000476837158 | 08510838+1147121
 516240 | Vilnius |           4873.0 | 3.5899999141693115 | 08513045+1148582
 516242 | Vilnius |           5061.0 | 3.7999999523162842 | 08513577+1153347
 516246 | Vilnius |           4946.0 | 3.5399999618530273 | 08514507+1147459
(8 rows)
"""


database.execute(
    """ UPDATE results
        SET passed_quality_control = false
        WHERE id in (504050, 504060, 504062, 504066, 516230, 516240, 516242, 516246)""")

# Has fixed it.

raise a
