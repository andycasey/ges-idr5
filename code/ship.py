
""" Fucking :shipit: """

import logging
import numpy as np
from astropy.io import fits
from astropy.table import Table
from datetime import datetime

import utils
from db import Database

logger = logging.getLogger("ges")


def wg_recommended_sp_template(database, input_path, output_path, wg,
    ext=-1, overwrite=False, **kwargs):
    """
    Produce a WG-recommended file of stellar parameters from the template 
    provided (on a per CNAME) basis.

    :param database:
        A database to connect to.

    :param input_path:
        The local file path of the WG-recommended template file.

    :param output_path:
        The local file path where to save the WG-recommended file to.

    :param wg:
        The working group.

    :param ext: [optional]
        The extension of the template to use for updating.

    :param overwrite: [optional]
        Overwrite the `output_path` if it already exists.
    """

    columns = [
        "teff",
        "e_teff",
        "e_pos_teff",
        "e_neg_teff",
        "nn_nodes_teff",
        "nn_spectra_teff",
        "enn_teff",
        "nne_teff",
        "sys_err_teff",
        
        "logg",
        "e_logg",
        "e_pos_logg",
        "e_neg_logg",
        "nn_nodes_logg",
        "nn_spectra_logg",
        "enn_logg",
        "nne_logg",
        "sys_err_logg",
        "lim_logg",
        
        "feh",
        "e_feh",
        "e_pos_feh",
        "e_neg_feh",
        "nn_nodes_feh",
        "nn_spectra_feh",
        "enn_feh",
        "nne_feh",
        "sys_err_feh",
        
        "xi",
        "e_xi",
        "e_pos_xi",
        "e_neg_xi",
        "nn_nodes_xi",
        "nn_spectra_xi",
        "enn_xi",
        "nne_xi",
        
        "mh",
        "e_mh",
        "e_pos_mh",
        "e_neg_mh",
        "nn_nodes_mh",
        "nn_spectra_mh",
        "enn_mh",
        "nne_mh",

        "alpha_fe",
        "e_alpha_fe",
        "e_pos_alpha_fe",
        "e_neg_alpha_fe",
        "nn_nodes_alpha_fe",
        "nn_spectra_alpha_fe",
        "enn_alpha_fe",
        "nne_alpha_fe",

        "tech",
        "remark",
        "peculi"
    ]

    translations = {
        "nn_nodes_teff": "NN_TEFF",
        "nn_nodes_logg": "NN_LOGG",
        "nn_nodes_feh": "NN_FEH",
        "nn_nodes_mh": "NN_MH",
        "nn_nodes_xi": "NN_XI",
        "nn_nodes_alpha_fe": "NN_ALPHA_FE"
    }
    default_translator = lambda x: x.upper()

    updated_data = {}
    for column in columns:
        updated_data[column] = []

    image = fits.open(input_path)

    if kwargs.pop("skip_verify", False):
        if image[0].header["NODE1"].strip() != "WG{}".format(wg):
            raise ValueError(
                "Template expected {} in NODE1 keyword, but WG{} given. "
                "Use skip_verify=False to override this error.".format(
                    image[0].header["NODE1"].strip(), wg))

    # For each cname in the template, fill in the entries from the database.
    N = len(image[ext].data)
    for i, cname in enumerate(image[ext].data["CNAME"]):

        logger.info("At row {}/{}: {}".format(i + 1, N, cname))

        record = database.retrieve_table(
            """ SELECT {columns}
                  FROM wg_recommended_results
                 WHERE wg = '{wg}'
                   AND cname = '{cname}'
            """.format(wg=wg, cname=cname, columns=", ".join(columns)))

        assert record is not None
        assert len(record) == 1

        for column in columns:
            updated_data[column].append(record[column][0])

    J = len(columns)
    for j, column in enumerate(columns):

        # Translate the column.
        fits_column = translations.get(column, None)
        if fits_column is None:
            fits_column = default_translator(column)

        logger.info("Updating column {}/{} in template: {} -> {}".format(
            j + 1, J, column, fits_column))

        if fits_column not in image[ext].data.dtype.names:
            logger.warn(
                "Column '{}' (from {}) not in template. Skipping..".format(
                    fits_column, column))
            continue
        
        try:
            image[ext].data[fits_column] = updated_data[column]
        
        except:
            logger.exception(
                "Could not update column '{}':".format(fits_column))
            continue

    # Update the release date.
    now = datetime.now()
    image[0].header["DATETAB"] = "{year}-{month}-{day}".format(
        year=now.year, month=now.month, day=now.day)
    image.writeto(output_path, clobber=overwrite)

    logger.info("Written WG{}-recommended file to {}".format(wg, output_path))

    return None