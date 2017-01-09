
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

        if record is None:
            # Keep whatever is currently there in the template.
            for column in columns:
                if default_translator(column) in image[ext].data.dtype.names:
                    updated_data[column].append(image[ext].data[default_translator(column)][i])
                else:
                    updated_data[column].append(np.nan)
        
        else:
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
       
        # Check for non-finite values and integer requirements.
        types = list(set(map(type, updated_data[column])))
        if len(types) == 2 \
        and any([isinstance(v, int) for v in updated_data[column]]) \
        and not np.all(np.isfinite(updated_data[column])):
            updated_data[column] = np.array(updated_data[column])
            bad = ~np.isfinite(updated_data[column])
            updated_data[column][bad] = 0
            updated_data[column] = updated_data[column].astype(int)
    
        try:
            image[ext].data[fits_column] = updated_data[column]
        
        except:
            logger.exception(
                "Could not update column '{}':".format(fits_column))

    # It's stupid that we should ever have to do this.
    propagate_columns = {
        "NN_TEFF": ("ENN_TEFF", "NNE_TEFF"),
        "NN_LOGG": ("ENN_LOGG", "NNE_LOGG"),
        "NN_FEH": ("ENN_FEH", "NNE_FEH"),
    }
    for column, propagate_to_columns in propagate_columns.items():
        for propagated_column in propagate_to_columns:
            image[ext].data[propagated_column] = image[ext].data[column]

    # If this is WG 11, then we will do the TECH flags.
    if wg == 11:
        logger.info("Doing TECH flags separately")
        max_length = -1
        concatenated_tech = []
        for i, cname in enumerate(image[ext].data["CNAME"]):

            logger.info("At row {}/{}: {}".format(i + 1, N, cname))

            if np.isfinite(updated_data["teff"][i]):
                record = database.retrieve_table(
                    """ SELECT w.id, string_agg(DISTINCT r.tech, '|') AS tech
                        FROM wg_recommended_results AS w,
                             results as r
                        WHERE w.wg = %s
                          AND r.id = ANY(
                                  w.provenance_ids_for_teff ||
                                  w.provenance_ids_for_logg ||
                                  w.provenance_ids_for_feh)
                          AND r.tech <> 'NaN'
                          AND r.tech <> ''
                          AND w.cname = %s
                        GROUP BY w.id;
                    """, (wg, cname))

            else:
                record = database.retrieve_table(
                    """ SELECT string_agg(DISTINCT r.tech, '|') AS tech
                        FROM results as r,
                             nodes as n
                        WHERE r.node_id = n.id
                          AND n.wg = %s
                          AND r.cname = %s
                          AND r.tech <> 'NaN'
                          AND r.tech <> '';
                    """, (wg, cname))

            if record is None:
                concatenated_tech.append("")

            else:
                tech = "|".join(sorted(map(str.strip, list(set(record["tech"][0].split("|"))))))
                concatenated_tech.append(tech)

                max_length = max([max_length, len(tech)])

            print("TECH", i, cname, concatenated_tech[-1])

        # Fuck this shit.
        if not "FixedLength" in input_path:
            cols = fits.ColDefs(
                [col for col in image[ext].columns[:-1]] \
                    + [fits.Column(name="TECH", format="A{}".format(max_length),
                                   array=concatenated_tech)])
            image[ext] = fits.BinTableHDU.from_columns(cols)

        else:
            if max_length > 250:
                logger.warn("Some TECH flags are going to be truncated!")
            image[ext].data["TECH"] = concatenated_tech


    # Update the release date.
    now = datetime.now()
    image[1].header["EXTNAME"] = "WGParametersWGAbundances"
    image[0].header["DATETAB"] = "{year}-{month}-{day}".format(
        year=now.year, month=now.month, day=now.day)

    # Create a temporary HDU
    hdu = fits.BinTableHDU()
    hdu.header.update({
        "RELEASE": image[0].header["RELEASE"],
        "DATETAB": image[0].header["DATETAB"],
        "INSTRUME": image[0].header["INSTRUME"],
        "NODE1": "WG{}".format(wg),
        "EXTNAME": "WGParametersWGAbundancesAdd",
        "EXTDUMMY": True,
        "COMMENT": "GES WG Recommended Parameters and Abundances"
    })
    if "SIMPLE" in hdu.header:
        del hdu.header["SIMPLE"]
    hdu.verify("fix")

    for column_name in ("ENN_TEFF", "ENN_LOGG", "ENN_FEH"):
        image[1].data[column_name] = image[1].data[column_name.replace("ENN_", "E_")]
        
    for column_name in ("NN_TEFF", "NN_LOGG", "NN_FEH", "NNE_TEFF", "NNE_LOGG", "NNE_FEH"):
        no_results = image[1].data[column_name] == 0
        image[1].data[column_name][no_results] = -1

    # Update with other nodes that contributed.
    contributed_nodes = database.retrieve_table(
        """ SELECT DISTINCT ON (n.name) n.name
            FROM nodes as n
            WHERE n.wg = %s 
              AND EXISTS(SELECT 1 FROM results AS r
                         WHERE r.node_id = n.id 
                           AND r.passed_quality_control)
            ORDER BY n.name ASC""", (wg, ))

    for i, node_name in enumerate(contributed_nodes["name"]):
        image[0].header["NODE{:.0f}".format(i + 2)] = node_name

    image.append(hdu)
    image.writeto(output_path, clobber=overwrite)
    
    logger.info("Written WG{}-recommended file to {}".format(wg, output_path))

    return None
