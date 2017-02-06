
"""
Update the benchmark parameters to include some values -- even if they are
uncertain -- and to include a less-biased value for HD 140283.
"""

from astropy.table import Table

input_path = "../fits-templates/benchmarks/GES_iDR5_FGKMCoolWarm_Benchmarks_AcceptedParams_01082016.fits"
output_path = "../fits-templates/benchmarks/GES_iDR5_FGKM_Benchmarks_ARC_29092016.fits"
overwrite = False


benchmarks = Table.read(input_path)
print("Read in benchmarks from {}".format(input_path))

updated_values = {
    "HD140283": {
        "TEFF": 5700,
        "E_TEFF": 200,
        "LOGG": 3.58, 
        "E_LOGG": 0.11,
        "FEH": -2.43,
    },
    "HD220009": {
        "TEFF": 4217,
        "E_TEFF": 60,
        "LOGG": 1.43,
        "E_LOGG": 0.12,
        "FEH": -0.75, 
    }
}

for ges_fld, params in updated_values.items():

    match = np.array([each.strip() == ges_fld for each in benchmarks["GES_FLD"]])

    for key, value in params.items():
        benchmarks[key][match] = value
        print("Updated {} = {} for {}".format(key, value, ges_fld))


benchmarks.write(output_path, overwrite=overwrite)
print("Written new file to {}".format(output_path))

