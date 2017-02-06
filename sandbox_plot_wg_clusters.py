
import yaml
from glob import glob

from code import (GESDatabase, plot, summary)
from astropy.table import Table

# Create a database object.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)

# Clean up bits and pieces...

cluster_names = ("M15", "NGC4372", "NGC4833", "M2", "NGC1904", "NGC6752",
    "M12", "NGC1261", "NGC362", "NGC1851", "NGC2808", "NGC104", "NGC5927",
    "NGC2243", "NGC6553", "NGC3532", "NGC6705", "NGC6253")


velocity_constraints = {
    "M15":      (-125, -90),  
    "NGC4372":  (60, 90),  
    "NGC4833":  (185, 215),  
    "M2":       (-20, 20),  
    "NGC1904":  (195, 215),  
    "NGC6752":  (-45, -10),  
    "M12":      (-50, -30),  
    "NGC1261":  (65, 85),  
    "NGC362 ":  (205, 240),  
    "NGC1851":  (305, 335),  
    "NGC2808":  (80, 130),  
    "NGC104":   (-33, 0),  
    "NGC5927":  (-115, -90),  
    "NGC2243":  (55, 65),  
    "NGC6553":  (-20, 20),  
    "NGC3532":  (-10, 10),  
    "NGC6705":  (-20, 50),  
    "NGC6253":  (-35, -20),  
}

metallicity_constraints = {
    "M15": (-5, -0.50),
    "NGC4372": (-3, -2.0),
    #"NGC4833": (-2.50, -1.75),
    "M2": (-5, -1.30), #--> check
    #-->"NGC1904": (-1.70, -1.00),
    "NGC6752": (-1.60, -1.30), # --> WTF is with the metal-poor stars here?
    "M12": (-1.60, -0.60),
    #"NGC1261": (-1.35, -0.80),
    "NGC362": (-5, -0.50),
    #"NGC1851": (-1.35, -0.70),
    "NGC2808": (-1.40, -0.60),
    "NGC104": (-0.85, -0.30),
    "NGC5927": (-0.55, -0.10 ),
    "NGC2243": (-1.00, -0.20),
    #"NGC6553": (-0.40, 0.10 ),
    "NGC3532": (-0.20, 0.30 ),
    "NGC6705": (-0.50, 0.50 ),
    "NGC6253": (-0.50 , 0.80),
}

"""
fig = plot.wg_metallicity_overview(database, x_axis="logg",
    cluster_names=cluster_names, velocity_constraints=velocity_constraints,
    metallicity_constraints=metallicity_constraints, wgs=(10, 11, 12, 13),
    ptp_ylim=None, sql_constraint="((s.wg = 10 and r.nn_nodes_teff > 1 and r.snr > 10 and r.teff > 4000 and r.teff < 6500) or (s.wg != 10))")
"""


isochrone_filenames = {
    "M15":     "isochrones/M15_Parsec_12Gyr_Z0.0001.dat",
    "NGC4372": "isochrones/NGC4372_Parsec2.9_10Gyr_Z0.0001.dat",
    "M2": "isochrones/M2_Parsec2.9_12Gyr_Z0.00034.dat",
    "NGC1904": "isochrones/NGC1904_Parsec2.9_12Gyr_Z0.00038.dat",
    "NGC6752": "isochrones/NGC6752_Parsec2.9_12Gyr_Z0.000438.dat",
    "M12": "isochrones/M12_Parsec2.9_12Gyr_Z0.000648.dat",
    "NGC1261": "isochrones/NGC1261_Parsec2.9_12Gyr_Z0.000816.dat",
    "NGC362": "isochrones/NGC362_Parsec2.9_12Gyr_Z0.000816.dat",
    "NGC1851": "isochrones/NGC1851_Parsec2.9_12Gyr_Z0.001.dat",
    "NGC2808": "isochrones/NGC2808_Parsec2.9_12Gyr_Z0.0011.dat",
    "NGC104": "isochrones/NGC104_Parsec2.9_12Gyr_Z0.002896.dat",
    "NGC5927": "isochrones/NGC5927_Parsec2.9_12Gyr_Z0.00491.dat",
    "NGC6553": "isochrones/NGC6553_Parsec2.9_12Gyr_Z0.01.dat",
    "NGC4833": "isochrones/NGC4833_Parsec_12Gyr_Z0.0002.dat",
    "NGC5927": "isochrones/NGC5927_Parsec_11Gyr_Z0.004.dat",
    "NGC2243": "isochrones/NGC2243_Parsec_4.5Gyr_Z0.004.dat",
    "NGC6705": "isochrones/NGC6705_Parsec_0.3Gyr_Z0.018.dat",
}  

fig = plot.hertzsprung_russell_diagrams(database, figsize=(10.725, 12.4125),
    cluster_names=cluster_names, velocity_constraints=velocity_constraints,
    metallicity_constraints=metallicity_constraints, wgs=(10, 11, ),
    isochrone_filenames=isochrone_filenames,
    sql_constraint="(r.wg = 10 and trim(r.setup) = 'HR10|HR21' and r.snr > 10 and r.nn_nodes_teff > 1) or (r.wg = 11 and r.nn_nodes_teff > 5)")#="((s.wg = 110 and r.nn_nodes_teff > 1 and r.snr > 10 and r.teff > 4000 and r.teff < 6500) or (s.wg != 110 and s.wg != 11) or (s.wg = 11 and r.nn_nodes_teff > 2))")
fig.savefig("figures/ges-dr5-clusters-all-wgs-hrds-20170206.pdf", dpi=300)
