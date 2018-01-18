# -*- coding: utf-8 -*-
""" This module contains tests enrichment of BEL graphs using Bio2BEL KEGG"""

from tests.constants import DatabaseMixin, enrichment_graph


class TestEnrich(DatabaseMixin):
    """Tests the enrichment of module"""

    def test_get_pathway_graph(self):
        graph = self.manager.get_pathway_graph('path:hsa00030')

        self.assertEqual(15, len(graph.number_of_nodes()))  # 14 proteins + pathway node
        self.assertEqual(14, len(graph.number_of_edges()))  # 14 edges protein -- pathway

    def test_enrich_kegg_pathway(self):
        graph_example = enrichment_graph()

        enriched_graph = self.manager.enrich_kegg_pathway(graph_example)

        self.assertEqual(16, len(
            enriched_graph.number_of_nodes()))  # 14 proteins in the pathway + gene of one of the proteins + pathway node
        self.assertEqual(15, len(enriched_graph.number_of_edges()))  # 14 edges protein -- pathway + gene -- pathway

    def test_enrich_kegg_protein(self):
        graph_example = enrichment_graph()

        enriched_graph = self.manager.enrich_kegg_protein(graph_example)

        self.assertEqual(5, len(enriched_graph.number_of_nodes()))  # 2 proteins + gene + pathway + new pathway
        self.assertEqual(4, len(enriched_graph.number_of_edges()))  # 3 edges + new one
