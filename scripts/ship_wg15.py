
"""
Ship a WG15 recommended SP product.
"""

import yaml
import logging

from code import GESDatabase, survey

# Initialize logging.
logger = logging.getLogger("ges")

# Create a database object.
db_filename = "db.yaml"
with open(db_filename, "r") as fp:
    credentials = yaml.load(fp)
database = GESDatabase(**credentials)


# Produce a WG15/Survey level homogenisation file
foo = survey.homogenise_survey_results(database)