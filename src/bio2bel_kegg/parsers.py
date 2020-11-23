# -*- coding: utf-8 -*-

"""This module parsers the KEGG pathway names file.

The "Complete list of pathways" file maps the KEGG identifiers to their corresponding pathway name .
"""

from typing import Optional

import pandas as pd

from bio2bel.utils import ensure_path
from .constants import KEGG_HUMAN_PATHWAYS_URL, KEGG_ORGANISM_URL, MODULE_NAME, PROTEIN_PATHWAY_HUMAN_URL

__all__ = [
    'get_pathway_df',
    'get_entity_pathway_df',
    'get_organisms_df',
    'parse_entry_line',
    'remove_first_word',
    'get_first_word',
    'parse_pathway_line',
    'parse_link_line',
    'parse_description',
    'get_description_properties',
    'kegg_properties_to_models',
    'process_protein_info_to_model'
]


def get_pathway_df(url: Optional[str] = None) -> pd.DataFrame:
    """Convert tab separated txt files to pandpathway = parse_pathway_lines(pathway_lines)as Dataframe.

    :param url: url from KEGG tab separated file
    :return: dataframe of the file
    """
    df = pd.read_csv(
        url or ensure_path(MODULE_NAME, KEGG_HUMAN_PATHWAYS_URL, path='pathways.tsv'),
        sep='\t',
        header=None,
        names=['kegg_pathway_id', 'name'],
    )
    # df['kegg_pathway_id'] = df['kegg_pathway_id'].map(_remove_path_prefix)
    return df


def get_entity_pathway_df(url: Optional[str] = None) -> pd.DataFrame:
    """Convert tab separated text files in to DataFrame.

    :param url: An optional url from a KEGG TSV file
    """
    df = pd.read_csv(
        url or ensure_path(MODULE_NAME, PROTEIN_PATHWAY_HUMAN_URL, path='protein_pathway.tsv'),
        sep='\t',
        header=None,
        names=['kegg_protein_id', 'kegg_pathway_id'],
    )
    # df['kegg_pathway_id'] = df['kegg_pathway_id'].map(_remove_path_prefix)
    return df


def get_organisms_df(url: Optional[str] = None) -> pd.DataFrame:
    """Convert tab separated txt files to pandas Dataframe.

    :param url: url from KEGG tab separated file
    :return: dataframe of the file
    :rtype: pandas.DataFrame
    """
    df = pd.read_csv(
        url or ensure_path(MODULE_NAME, KEGG_ORGANISM_URL, path='organisms.tsv'),
        sep='\t',
        header=None,
        names=[
            'kegg_id',
            'kegg_code',
            'name',
            # fourth column is the taxonomy hierarchy
        ],
        usecols=[0, 1, 2],
    )
    df['common_name'] = df['name'].map(lambda name: name.replace(')', '').split(' (')[1].capitalize() if len(name.replace(')', '').split(' (')) > 1 else '')
    df['name'] = df['name'].map(lambda name: name.replace(')', '').split(' (')[0].capitalize())
    return df


# -*- coding: utf-8 -*-

"""This module parsers the description files -> http://rest.kegg.jp/get/ in KEGG RESTful API."""

import re

from requests import Response

from bio2bel_kegg.constants import DBLINKS, PROTEIN_RESOURCES


def parse_entry_line(line):
    """Parse entry line to tuple.

    :param line:
    :rtype tuple
    :return: tuple of entry
    """
    return tuple(
        line.strip(' ')
        for line in line.split()[1:]
    )


def remove_first_word(string):
    """Remove the first word of the line.

    :param str string: string
    :rtype str
    :return: string without the first word
    """
    return string.split(' ', 1)[1].strip()


def get_first_word(string):
    """Get the first word of the line.

    :param str string: string
    :rtype str
    :return: string with the first word
    """
    return string.split(' ', 1)[0]


def parse_pathway_line(line):
    """Parse entry pathway line to tuple.

    :param line:
    :rtype tuple
    :return: tuple of entry
    """
    line = remove_first_word(line)

    return tuple(
        line.strip(' ')
        for line in re.split(r'\s{2,}', line)
    )


def parse_link_line(line):
    """Parse entry dblink line to tuple.

    :param line:
    :rtype tuple
    :return: tuple of entry
    """
    line = remove_first_word(line)

    column, link_id = line.split(":")

    return column.strip(), link_id.strip()


def parse_description(response: Response):
    """Parse the several properties in the description file given an KEGG identifier using the KEGG API.

    :rtype: dict
    :return: description dictionary
    """
    description = {}

    for line in response.iter_lines():
        line = line.decode('utf-8')

        if not line.startswith(' '):
            keyword = get_first_word(line)

        if keyword == 'ENTRY':
            description['ENTRY'] = parse_entry_line(line)

        elif keyword == 'NAME':
            entry_name = parse_entry_line(line)
            if entry_name:
                # If there is a name, take the first element of the tuple and strip semi colon
                # in case there are multiple names
                description['ENTRY_NAME'] = entry_name[0].strip(';')

        elif keyword == 'PATHWAY':

            if 'PATHWAY' not in description:
                description['PATHWAY'] = [parse_pathway_line(line)]
            else:
                description['PATHWAY'].append(parse_pathway_line(line))

        elif keyword == 'DBLINKS':

            if 'DBLINKS' not in description:
                description['DBLINKS'] = [parse_link_line(line)]
            else:
                description['DBLINKS'].append(parse_link_line(line))

    return description


def get_description_properties(description, description_property, columns):
    """Get specific description properties.

    :param dict protein_description: id for the query
    :param str description_property: main property in the description
    :param list columns: columns to be filtered
    :rtype: dict
    :return: description dictionary
    """
    return {
        pair[0]: pair[1]
        for pair in description[description_property]
        if pair[0] in columns
    }


def kegg_properties_to_models(kegg_attributes):
    """Modify the kegg attribute dictionary to match the db '{}_id' formatting.

    :param dict kegg_attributes: kegg description dictionary
    :rtype: dict
    :return: dictionary with bio2bel_kegg adapted keys
    """
    return {
        '{}_id'.format(key.lower()): value
        for key, value in kegg_attributes.items()
        if len(value) < 255
    }


def process_protein_info_to_model(response: Response):
    """Process description.

    :param response: response from KEGG API
    :type: dict
    :return: protein model attributes
    """
    # Get protein description from KEGG API
    description = parse_description(response)
    # Filters out db link columns
    protein_as_dict = get_description_properties(
        description=description,
        description_property=DBLINKS,
        columns=PROTEIN_RESOURCES
    )
    # Adapt the dict keys to match protein model columns
    return kegg_properties_to_models(protein_as_dict)
