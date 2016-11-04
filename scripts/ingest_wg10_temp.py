
"""
Import WG recommended results on a per-SETUP basis as if they were a per-CNAME
basis.

NOTE: THIS IS A HACK AND IS ONLY HERE TO TEMPORARILY BE ABLE TO CREATE A WG15
      FILE!
"""

import logging
import numpy as np
import psycopg2 as pg
import yaml
from glob import glob
from astropy.table import Table

from code import GESDatabase

logger = logging.getLogger("ges")

# Get credentials.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)

# Create a database object.
database = GESDatabase(**credentials)

# Ingest shitty WG10 results
data = Table.read("recommended-results/GES_iDR5_WG10_Recommended_PERSETUP.fits")

default_row = { "wg": 10 }

columns = ("wg", # For default row
    "cname", "filename", "setup", "snr",
    "vel", "e_vel", "vrot", "e_vrot",     
    "teff", "e_teff", "nn_teff", "enn_teff", "nne_teff", "sys_err_teff",
    "logg", "e_logg", "nn_logg", "enn_logg", "nne_logg", "sys_err_logg", "lim_logg",
    "feh", "e_feh", "nn_feh", "enn_feh", "nne_feh", "sys_err_feh",
    "xi", "e_xi", "nn_xi", "enn_xi", "nne_xi",
    "mh", "e_mh", "nn_mh", "enn_mh", "nne_mh",
    "alpha_fe", "e_alpha_fe", "nn_alpha_fe", "enn_alpha_fe", "nne_alpha_fe", 
    "vrad", "e_vrad", "vsini", "e_vsini",
    "peculi", "remark", "tech",

    "lim_vsini",
    "teff_phot",
    "e_teff_phot",
    "teff_irfm",
    "e_teff_irfm",
    "fbol_irfm",

    "spt",
    "veil",
    "e_veil",
    "ew_li",
    "lim_ew_li",
    "e_ew_li",
    "ewc_li",
    "lim_ewc_li",
    "e_ewc_li",
    "ew_ha_acc",
    "e_ew_ha_acc",
    "ha10",
    "e_ha10",
    "ew_ha_chr",
    "e_ew_ha_chr",
    "fha_chr",
    "e_fha_chr",
    "fwzi",
    "e_fwzi",
    "ew_hb_chr",
    "e_ew_hb_chr",
    "fhb_chr",
    "e_fhb_chr",
    "log_mdot_acc",
    "e_log_mdot_acc",
    "log_l_acc",
    "e_log_l_acc",
    "gamma",
    "e_gamma",
    "convol",
    "e_convol",
    "m_alpha",
    "m_grid",
    "m_broad",
    "m_loops",
    "m_name",
)

# Update formats, as necessary.
tmp_key_format = "{}_NEW_DTYPE"
for key, new_dtype in _FITS_FORMAT_ADAPTERS.items():
    data[tmp_key_format.format(key.upper())] = np.array(data[key.upper()], dtype=new_dtype)
    del data[key.upper()]
    data.rename_column(tmp_key_format.format(key.upper()), key.upper())


data = data.group_by("CNAME")

for i, group in enumerate(data.groups):

    # Fuck it: take the first setup result for this star
    row = group[0]

    logger.info("Ingesting row {}/{} from WG{}".format(i, N, wg))
    row_data = {}
    row_data.update(default_row)
    row_data.update(dict(zip(columns[1:], [row[c.upper()] for c in columns[1:]])))

    # Trim strings!
    for k in row_data.keys():
        if isinstance(row_data[k], str):
            row_data[k] = row_data[k].strip()

    database.execute(
        "INSERT INTO wg_recommended_results ({}) VALUES ({})".format(
            ", ".join(columns),
            ", ".join(["%({})s".format(column) for column in columns])),
        row_data)

database.connection.commit()
