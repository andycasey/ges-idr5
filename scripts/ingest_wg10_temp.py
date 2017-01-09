
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



_FITS_FORMAT_ADAPTERS = {
    "snr": float,
    "vel": float,
    "e_vel": float,
    "vrot": float,
    "e_vrot": float,
    "teff": float,
    "e_teff": float,
    "nn_teff": int,
    "enn_teff": float,
    "nne_teff": float,
    "sys_err_teff": float,
    "logg": float,
    "e_logg": float,
    "nn_logg": int,
    "enn_logg": float,
    "nne_logg": float,
    "sys_err_logg": float,
    "lim_logg": int,
    "feh": float,
    "e_feh": float,
    "nn_feh": int,
    "enn_feh": float,
    "nne_feh": float,
    "sys_err_feh": float,
    "xi": float,
    "e_xi": float,
    "nn_xi": int,
    "enn_xi": float,
    "nne_xi": float,
    "mh": float,
    "e_mh": float,
    "nn_mh": int,
    "enn_mh": float,
    "nne_mh": float,
    "alpha_fe": float,
    "e_alpha_fe": float,
    "nn_alpha_fe": int,
    "enn_alpha_fe": float,
    "nne_alpha_fe": float,
    "vrad": float,
    "e_vrad": float,
    "vsini": float,
    "e_vsini": float,

    "lim_vsini": int,
    "teff_phot": float,
    "e_teff_phot": float,
    "teff_irfm": float,
    "e_teff_irfm": float,
    "fbol_irfm": float,

    "veil": float,
    "e_veil": float,
    "ew_li":  float,
    "lim_ew_li": int,
    "e_ew_li": float, 
    "ewc_li": float, 
    "lim_ewc_li": int, 
    "e_ewc_li": float,
    "ew_ha_acc": float,
    "e_ew_ha_acc": float,
    "ha10": float,
    "e_ha10": float,
    "ew_ha_chr": float,
    "e_ew_ha_chr": float,
    "fha_chr": float,
    "e_fha_chr": float,
    "fwzi": float,
    "e_fwzi": float,
    "ew_hb_chr": float,
    "e_ew_hb_chr": float,
    "fhb_chr": float,
    "e_fhb_chr": float,
    "log_mdot_acc": float,
    "e_log_mdot_acc": float,
    "log_l_acc": float,
    "e_log_l_acc": float,
    "gamma": float,
    "e_gamma": float,
    "convol": float,
    "e_convol": float,
    "m_alpha": float,
    "m_broad": float,
    "m_loops": float,
}

wg = 110
default_row = { "wg": wg }

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


data = data.group_by(["CNAME", "SETUP"])
N = len(data.groups)

for i, group in enumerate(data.groups):

    row = group[0]
    """
    hr1021 = np.where(group["SETUP"] == "HR10|HR21")[0]
    if not hr1021:
        row = group[0]
    else:
        row = group[hr1021[0]]
    """

    row_data = {}
    row_data.update(default_row)
    row_data.update(dict(zip(columns[1:], [row[c.upper()] for c in columns[1:]])))

    logger.info("Ingesting row {}/{} from WG{}: {} / {} / {}".format(i, N, row_data["wg"], 
        row_data["teff"], row_data["logg"], row_data["feh"]))

  
    # Trim strings!
    for k in row_data.keys():
        if isinstance(row_data[k], str):
            row_data[k] = row_data[k].strip()

    database.execute(
        "INSERT INTO wg_recommended_results ({}) VALUES ({})".format(
            ", ".join(columns),
            ", ".join(["%({})s".format(column) for column in columns])),
        row_data)

    if (i % 100) == 0:
        print("COMMITTING")
        database.connection.commit()

database.connection.commit()
