"""Unit tests for the reusable CGR implementation."""

import tempfile
import unittest
from pathlib import Path

import numpy as np

from cgr import (
    count_kmers,
    fcgr_distance,
    fcgr_matrix,
    kmer_probabilities,
    read_fasta,
)


class TestKmerUtilities(unittest.TestCase):
    def test_counts_overlapping_kmers(self):
        self.assertEqual(count_kmers("AAAA", 2), {"AA": 3})

    def test_ambiguous_windows_are_skipped(self):
        self.assertEqual(count_kmers("AANAA", 2), {"AA": 2})

    def test_k_must_be_positive(self):
        with self.assertRaises(ValueError):
            count_kmers("ACGT", 0)

    def test_probabilities_use_occurrence_total(self):
        probabilities = kmer_probabilities({"AA": 3, "AC": 1})
        self.assertAlmostEqual(probabilities["AA"], 0.75)
        self.assertAlmostEqual(probabilities["AC"], 0.25)
        self.assertAlmostEqual(sum(probabilities.values()), 1.0)

    def test_empty_counts_return_empty_probabilities(self):
        self.assertEqual(kmer_probabilities({}), {})


class TestFCGR(unittest.TestCase):
    def test_single_base_corner_assignment(self):
        matrix = fcgr_matrix({"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}, 1)
        expected = np.array([[0.25, 0.25], [0.25, 0.25]])
        np.testing.assert_allclose(matrix, expected)

    def test_matrix_shape_and_probability_mass(self):
        probabilities = kmer_probabilities(count_kmers("ACGTACGT", 2))
        matrix = fcgr_matrix(probabilities, 2)
        self.assertEqual(matrix.shape, (4, 4))
        self.assertAlmostEqual(float(matrix.sum()), 1.0)

    def test_distance(self):
        first = np.zeros((2, 2))
        second = np.ones((2, 2))
        self.assertAlmostEqual(fcgr_distance(first, second), 2.0)
        self.assertEqual(fcgr_distance(first, first), 0.0)

    def test_distance_rejects_different_shapes(self):
        with self.assertRaises(ValueError):
            fcgr_distance(np.zeros((2, 2)), np.zeros((4, 4)))

    def test_fasta_reader_skips_headers_and_blank_lines(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta_path = Path(directory) / "sample.fasta"
            fasta_path.write_text(">sample\nACGT\n\nTGCA\n", encoding="utf-8")
            self.assertEqual(read_fasta(fasta_path), "ACGTTGCA")


if __name__ == "__main__":
    unittest.main()
