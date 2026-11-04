"""Core DNA sequence and Frequency Chaos Game Representation utilities."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Mapping, Union

import numpy as np

PathLike = Union[str, Path]

# Bit coordinates for the conventional CGR nucleotide corner assignment.
_NUCLEOTIDE_BITS = {
    "A": (0, 0),
    "C": (0, 1),
    "G": (1, 1),
    "T": (1, 0),
}


def read_fasta(path: PathLike) -> str:
    """Read and concatenate sequence lines from a FASTA file.

    Header lines (those beginning with '>') and blank lines are ignored.
    Sequence characters are returned in uppercase; validation is performed
    when k-mers are counted.
    """
    sequence_parts = []
    with Path(path).open("r", encoding="utf-8") as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line and not line.startswith(">"):
                sequence_parts.append(line)
    return "".join(sequence_parts).upper()


def count_kmers(sequence: str, k: int) -> Counter[str]:
    """Count valid, overlapping DNA k-mers in *sequence*.

    Windows containing ambiguous or non-DNA characters are skipped.
    """
    _validate_k(k)
    normalized = sequence.upper()
    counts: Counter[str] = Counter()

    for start in range(max(0, len(normalized) - k + 1)):
        kmer = normalized[start : start + k]
        if all(base in _NUCLEOTIDE_BITS for base in kmer):
            counts[kmer] += 1

    return counts


def kmer_probabilities(kmer_counts: Mapping[str, int]) -> dict[str, float]:
    """Normalize k-mer counts by the total number of valid occurrences."""
    if any(not isinstance(count, (int, np.integer)) or count < 0
           for count in kmer_counts.values()):
        raise ValueError("k-mer counts must be non-negative integers")

    total = sum(kmer_counts.values())
    if total == 0:
        return {}

    return {kmer: count / total for kmer, count in kmer_counts.items() if count}


def fcgr_matrix(probabilities: Mapping[str, float], k: int) -> np.ndarray:
    """Build a 2**k by 2**k FCGR matrix from k-mer probabilities.

    Matrix rows correspond to the CGR y-coordinate and columns to x.
    Multiple contributions to a cell are accumulated.
    """
    _validate_k(k)
    matrix_size = 2**k
    matrix = np.zeros((matrix_size, matrix_size), dtype=float)

    for kmer, probability in probabilities.items():
        if len(kmer) != k:
            raise ValueError(f"Expected k-mers of length {k}, got {kmer!r}")
        if not np.isfinite(probability) or probability < 0:
            raise ValueError("Probabilities must be finite and non-negative")

        x = 0
        y = 0
        for base in kmer.upper():
            if base not in _NUCLEOTIDE_BITS:
                raise ValueError(f"Unsupported nucleotide {base!r} in {kmer!r}")
            x_bit, y_bit = _NUCLEOTIDE_BITS[base]
            x = (x << 1) | x_bit
            y = (y << 1) | y_bit

        matrix[y, x] += probability

    return matrix


def fcgr_distance(first: np.ndarray, second: np.ndarray) -> float:
    """Return the Euclidean distance between two same-shaped FCGR matrices."""
    first_array = np.asarray(first, dtype=float)
    second_array = np.asarray(second, dtype=float)

    if first_array.shape != second_array.shape:
        raise ValueError("FCGR matrices must have identical shapes")
    if not np.isfinite(first_array).all() or not np.isfinite(second_array).all():
        raise ValueError("FCGR matrices must contain only finite values")

    return float(np.linalg.norm(first_array - second_array))


def _validate_k(k: int) -> None:
    if isinstance(k, bool) or not isinstance(k, int) or k < 1:
        raise ValueError("k must be a positive integer")
