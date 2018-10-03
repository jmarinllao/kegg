# -*- coding: utf-8 -*-

"""This module contains multiple parsers for the KEGG public data sources"""

from . import entities, pathways, description
from .description import *
from .entities import *
from .pathways import *

__all__ = (
        description.__all__ +
        entities.__all__ +
        pathways.__all__
)
