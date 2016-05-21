
""" A convenience object for databases. """

import logging
import psycopg2 as pg
from time import time

logger = logging.getLogger("ges")

class Database(object):

    def __init__(self, **kwargs):
        self.connection = pg.connect(**kwargs)
        return None


    def update(self, query, values=None, full_output=False, **kwargs):
        """
        Update the database with a SQL query.

        :param query:
            The SQL query to execute.

        :type query:
            str

        :param values: [optional]
            Values to use when formatting the SQL string.

        :type values:
            tuple or dict
        """

        logger.debug("Running SQL update query: {}".format(query))
        names, results, cursor = self.execute(query, values, **kwargs)
        return (names, results, cursor) if full_output else cursor.rowcount
        

    def retrieve(self, query, values=None, full_output=False, **kwargs):
        """
        Retrieve some data from the database.

        :param query:
            The SQL query to execute.

        :type query:
            str

        :param values: [optional]
            Values to use when formatting the SQL string.

        :type values:
            tuple or dict
        """

        names, results, cursor = self.execute(query, values, fetch=True,
            **kwargs)
        return (names, results, cursor.rowcount) if full_output else results


    def execute(self, query, values=None, fetch=False, **kwargs):
        """
        Execute some SQL from the database.

        :param query:
            The SQL query to execute.

        :type query:
            str

        :param values: [optional]
            Values to use when formatting the SQL string.

        :type values:
            tuple or dict
        """

        t_init = time()
        try:
            with self.connection.cursor() as cursor:
                cursor.execute(query, values)
                if fetch: results = cursor.fetchall()
                else: results = None

        except pg.ProgrammingError:
            logger.exception("SQL query failed: {0}, {1}".format(query, values))
            cursor.close()
            raise
        
        else:
            taken = 1e3 * (time() - t_init)
            try:
                logger.debug("Took {0:.0f} ms for SQL query {1}".format(taken,
                    " ".join((query % values).split())))
            except (TypeError, ValueError):
                logger.debug("Took {0:.0f} ms for SQL query {1} with values {2}"\
                    .format(taken, query, values))
        
        names = None if cursor.description is None \
            else tuple([column[0] for column in cursor.description])
        return (names, results, cursor)


    def node_id(self, description):
        """
        Return the identifer for a node based on its description. 

        If the `description` given is an integer, no search will occur and that
        value will be returned. If a string-like object is provided in the
        `description`, then the expected format is 'wg.node_query', or just
        'node_query' if the node only applies to one working group.

        :param description:
            The search term to use to identify the node.
        """

        try:
            node = int(description)
        except:
            # If '.' is in the node descriptor, 
            # Split the node descriptor by '.'

            if "." in description:
                wg, node_description = description.split(".")
            else:
                raise NotImplementedError

            # Search by the name, filter by wg?
            # TODO

        else:
            return node



