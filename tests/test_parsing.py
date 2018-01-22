# -*- coding: utf-8 -*-
""" This module contains tests for parsing KEGG files"""

from unittest import TestCase

from bio2bel_kegg.constants import DBLINKS, PROTEIN_RESOURCES
from bio2bel_kegg.models import Pathway, Protein
from bio2bel_kegg.parsers.description import parse_description, get_description_properties
from tests.constants import DatabaseMixin


class TestParse(DatabaseMixin):
    """Tests the parsing module"""

    def test_pathway_count(self):
        pathway_number = self.manager.session.query(Pathway).count()
        self.assertEqual(9, pathway_number)

    def test_protein_count(self):
        protein_number = self.manager.session.query(Protein).count()
        self.assertEqual(29, protein_number)

    def test_pathway_protein_1(self):
        pathway = self.manager.get_pathway_by_id('path:hsa00030')
        self.assertIsNotNone(pathway)
        self.assertEqual(14, len(pathway.proteins))

    def test_pathway_protein_2(self):
        pathway = self.manager.get_pathway_by_id('path:hsa00010')
        self.assertIsNotNone(pathway)
        self.assertEqual(16, len(pathway.proteins))

    def test_protein_pathway_1(self):
        protein = self.manager.session.query(Protein).filter(Protein.kegg_id == 'hsa:5214').one_or_none()
        self.assertIsNotNone(protein)
        self.assertEqual(2, len(protein.pathways))
        self.assertEqual(
            {'path:hsa00030', 'path:hsa00010'}, {
                pathway.kegg_id
                for pathway in protein.pathways
            }
        )


class TestDescriptionParse(TestCase):
    """ Description parsing test"""

    def test_description_protein(self):
        # Dictionary of list of tuples

        PFKP_protein = parse_description(identifier='hsa:5214')

        self.assertEqual(
            [('hsa00010', 'Glycolysis / Gluconeogenesis'),
             ('hsa00030', 'Pentose phosphate pathway'),
             ('hsa00051', 'Fructose and mannose metabolism'),
             ('hsa00052', 'Galactose metabolism'),
             ('hsa01100', 'Metabolic pathways'),
             ('hsa01200', 'Carbon metabolism'),
             ('hsa01230', 'Biosynthesis of amino acids'),
             ('hsa03018', 'RNA degradation'),
             ('hsa04152', 'AMPK signaling pathway'),
             ('hsa04919', 'Thyroid hormone signaling pathway'),
             ('hsa05230', 'Central carbon metabolism in cancer')
             ],
            PFKP_protein['PATHWAY']
        )

        self.assertEqual(
            [('NCBI-GeneID', '5214'),
             ('NCBI-ProteinID', 'NP_002618'),
             ('OMIM', '171840'),
             ('HGNC', '8878'),
             ('Ensembl', 'ENSG00000067057'),
             ('Vega', 'OTTHUMG00000017556'),
             ('Pharos', 'Q01813(Tbio)'),
             ('UniProt', 'Q01813'),
             ],
            PFKP_protein[DBLINKS]
        )

        description_links = get_description_properties(PFKP_protein, DBLINKS, PROTEIN_RESOURCES)

        self.assertDictEqual(
            {'HGNC': '8878', 'UniProt': 'Q01813'},
            description_links
        )