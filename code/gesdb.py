
""" A specialized database class for Gaia-ESO Survey data releases. """

import logging
import numpy as np
from astropy.io import fits
from astropy.table import Table

import utils
from db import Database

logger = logging.getLogger("ges")


class GESDatabase(Database):

    def __init__(self, *args, **kwargs):
        super(GESDatabase, self).__init__(*args, **kwargs)


    def create_or_retrieve_node_id(self, wg, node_name):
        """
        Reteive a unique identifier for a node, or create one.

        :param wg:
            The working group (e.g., 10).

        :param node_name:
            The name of the node.

        """

        try:
            return self.retrieve_node_id(wg, node_name)
        except UnknownNodeError:
            return self._create_node(wg, node_name)


    def retrieve_node_id(self, wg, node_name):
        """
        Retrieve a unique identifier for a node.

        :param wg:
            The working group (w.g., 10).

        :param node_name:
            The name of the node.

        :raises UnknownNodeError:
            If no node exists.

        :returns:
            The identifier.
        """

        result = self.retrieve("""SELECT id FROM nodes
            WHERE wg = %s AND lower(name) = %s""",
            (utils.wg_as_int(wg), node_name.strip().lower(), ))

        if not result:
            raise UnknownNodeError("node does not exist")
        else:
            return int(result[0][0])


    def _create_node(self, wg, node_name):
        """
        Create a node.

        :param wg:
            The working group (e.g., 10).

        :param node_name:
            The name of the node.
        """

        wg = utils.wg_as_int(wg)
        node_name = node_name.strip()

        result = self.execute(
            """INSERT INTO nodes (wg, name) VALUES (%s, %s) RETURNING id""",
            (wg, node_name), fetch=True)
        node_id = int(result[1][0][0])

        logger.info("Created node '{}' in WG{} with id {}".format(
            node_name, wg, node_id))
        return node_id


    def ingest_recommended_results_from_previous_dr(self, filename, extension=-1):
        """
        Ingest recommended results from a node FITS file.

        :param filename:
            A node template file in FITS format.

        :param extension: [optional]
            The extension index to read from.

        :returns:
            The number of rows inserted.
        """

        image = fits.open(filename)
        data = image[extension].data

        columns = ("cname", "ges_fld", "object", "filename", "ges_type",
            "teff", "e_teff", "logg", "e_logg", "mh", "e_mh", "xi", "e_xi",
            "peculi", "remark", "tech")

        fits_format_adapters = {
            "teff": float,
            "e_teff": float,
            "logg": float,
            "e_logg": float,
            "mh": float,
            "e_mh": float,
            "xi": float,
            "e_xi": float,
        }

        N = len(data)
        for i, row in enumerate(data):
            logger.info("Ingesting recommended row {}/{}".format(i, N))
            row_data = {}
            for column in columns:
                value = row[column]
                f = fits_format_adapters.get(column, None)
                if f is not None:
                    value = f(value)
                row_data[column] = value

            self.execute(
                "INSERT INTO recommended_idr4 ({}) VALUES ({})".format(
                    ", ".join(columns),
                    ", ".join(["%({})s".format(column) for column in columns])),
                row_data)

        self.connection.commit()
        return N


    def ingest_node_results(self, filename, extension=-1):
        """
        Ingest results from a node FITS file.

        :param filename:
            A node template file in FITS format.

        :param extension: [optional]
            The extension index to read from.

        :returns:
            The number of rows inserted.
        """

        # Which node is this?
        wg, node_name = utils.parse_node_filename(filename)
        node_id = self.retrieve_node_id(wg, node_name)

        # Start ingesting results.
        data = Table.read(filename, hdu=extension)

        default_row = {"node_id": node_id}
        columns = (
            "node_id", "cname", "filename", "setup", "snr",
            "vel", "e_vel", "vrot", "e_vrot",     
            "teff", "e_teff", "nn_teff", "enn_teff", "nne_teff", "sys_err_teff",
            "logg", "e_logg", "nn_logg", "enn_logg", "nne_logg", "sys_err_logg", "lim_logg",
            "feh", "e_feh", "nn_feh", "enn_feh", "nne_feh", "sys_err_feh",
            "xi", "e_xi", "nn_xi", "enn_xi", "nne_xi",
            "mh", "e_mh", "nn_mh", "enn_mh", "nne_mh",
            "alpha_fe", "e_alpha_fe", "nn_alpha_fe", "enn_alpha_fe", "nne_alpha_fe", 
            "vrad", "e_vrad", "vsini", "e_vsini",
            "peculi", "remark", "tech")

        fits_format_adapters = {
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
        }

        # Update formats, as necessary.
        tmp_key_format = "{}_NEW_DTYPE"
        for key, new_dtype in fits_format_adapters.items():
            data[tmp_key_format.format(key.upper())] = np.array(data[key.upper()], dtype=new_dtype)
            del data[key.upper()]
            data.rename_column(tmp_key_format.format(key.upper()), key.upper())

        N = len(data)
        for i, row in enumerate(data):
            logger.info("Ingesting row {}/{} from node WG{}: {}".format(i, N,
                wg, node_name))
            row_data = {}
            row_data.update(default_row)
            row_data.update(dict(zip(columns[1:], [row[c.upper()] for c in columns[1:]])))

            self.execute(
                "INSERT INTO results ({}) VALUES ({})".format(
                    ", ".join(columns),
                    ", ".join(["%({})s".format(column) for column in columns])),
                row_data)

        self.connection.commit()
        return N

        
    def ingest_spectra_masterlist(self, filename, extension=-1):
        """
        Ingest a master list of spectra from a FITS template file.

        :param filename:
            A FITS template file that contains the masterlist of all spectra.

        :returns:
            The number of rows inserted.
        """

        image = fits.open(filename)
        data = image[extension].data

        # Create mapper between FITS and database columns.
        columns = ("cname", "ges_fld", "object", "filename", "ges_type", "setup",
            "wg", "instrument", "ra", "dec", "snr", "vel", "e_vel", "vrot", 
            "e_vrot", "teff_irfm", "e_teff_irfm", "peculi", "remark", "tech")
        fits_column_adapters = {
            "instrument": "instrume"
        }
        fits_format_adapters = {
            "wg": utils.safe_int,
            "ra": float,
            "dec": float,
            "snr": float,
            "vel": float,
            "e_vel": float,
            "vrot": float,
            "e_vrot": float,
            "teff_irfm": float,
            "e_teff_irfm": float,
        }
        
        N = len(data)
        for i, row in enumerate(data):
            logger.info("Inserting row {}/{}".format(i, N))

            values = []
            for col in columns:
                use_col = fits_column_adapters.get(col, col)
                value = row[use_col]

                # Formatting.
                if col in fits_format_adapters:
                    f = fits_format_adapters[col]
                    value = f(value)
                values.append(value)

            self.execute(
                "INSERT INTO spectra ({}) VALUES ({})".format(
                    ", ".join(columns), ", ".join(["%s"] * len(columns))),
                values)

        self.connection.commit()

        return N



    def ingest_magrini_photometric_temperatures(self, filename, extension=-1):
        """
        Ingest a FITS table containing CNAMEs and photometric temperatures.

        :param filename:
            A FITS table.

        :param extension: [optional]
            The HDU extension that contains the photometric temperatures.
        """

        image = fits.open(filename)
        data = image[extension].data

        # The columns might be different, but in general if we lowerize them all
        # then we are looking for: 
        # ('CNAME_2', 'GES_FLD', 'teffjk', 'jk', 'FILENAME')
        cname_col, teff_col = (data.dtype.names[0], "teffjk")

        # Update the value in the spectra table, unless it already exists.    
        N = 0
        for row in data:
            result = self.execute(
                """ UPDATE spectra
                    SET teff_irfm = %s 
                    WHERE   cname = %s AND
                            teff_irfm = 'NaN'""",
                        (float(row[teff_col]), row[cname_col], ))
        return True


class UnknownNodeError(BaseException):
    pass
