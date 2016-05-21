
""" A specialized database class for Gaia-ESO Survey data releases. """

import logging
from astropy.io import fits

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
        image = fits.open(filename)
        data = image[extension].data

        default_row = {"node_id": node_id}
        columns = ("node_id", "cname", "filename", "setup",
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
            logger.info("Ingesting row {}/{} from node WG{}: {}".format(i, N,
                wg, node_name))
            row_data = {}
            row_data.update(default_row)
            for column in columns[1:]:
                value = row[column]
                f = fits_format_adapters.get(column, None)
                if f is not None:
                    value = f(value)
                row_data[column] = value

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


class UnknownNodeError(BaseException):
    pass
