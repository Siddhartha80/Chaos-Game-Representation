# Chaos Game Representation for DNA Sequences

A small research project exploring **Chaos Game Representation (CGR)** and **Frequency Chaos Game Representation (FCGR)** for visualizing DNA sequences and comparing sequence composition.

The repository includes the original research notebooks, example FASTA sequences, and generated CGR images. A lightweight Python package is also provided so the core sequence-processing operations can be reused and tested independently of the notebooks.

## Repository contents

| Path | Description |
| --- | --- |
| `cgr/` | Reusable Python implementation of FASTA parsing, k-mer counting, FCGR generation, and matrix distance |
| `tests/` | Unit tests for the reusable implementation |
| `CGR_SEQ_SIMILARITY.ipynb` | Notebook exploring k-mer FCGRs, sequence distances, and hierarchical clustering |
| `Leigh Syndrome.ipynb` | Notebook demonstrating k-mer-based FCGR analysis on a mitochondrial reference sequence |
| `DNA_SEQUENCES/` | Example DNA FASTA files and associated generated images |
| `ANIMAL_GENOME/` | Example animal genome FASTA files and associated generated images |
| `NC_012920.fasta` | FASTA sequence used by the Leigh Syndrome notebook |
| `BIO_Report_IBS3_RVM.pdf`, `BIO_IBS3_E4_RVM_12_JAN_2023.pptx` | Supporting project report and presentation |

## What is CGR?

Chaos Game Representation maps a nucleotide sequence to points in a unit square. Each nucleotide moves the current point halfway toward its corresponding corner:

| Nucleotide | Coordinate |
| --- | --- |
| A | (0, 0) |
| C | (0, 1) |
| G | (1, 1) |
| T | (1, 0) |

FCGR divides this representation into a grid and records k-mer frequencies in the corresponding cells. This provides a compact representation that can be visualized or compared across sequences.

## Getting started

Python 3.9 or later is recommended.

1. Clone the repository and enter its directory:

   ```bash
   git clone https://github.com/Siddhartha80/Chaos-Game-Representation.git
   cd Chaos-Game-Representation
   ```

2. Create and activate a virtual environment:

   ```bash
   python -m venv .venv
   ```

   On macOS/Linux:

   ```bash
   source .venv/bin/activate
   ```

   On Windows PowerShell:

   ```powershell
   .venv\Scripts\Activate.ps1
   ```

3. Install the project dependencies:

   ```bash
   python -m pip install --upgrade pip
   pip install -r requirements.txt
   ```

4. Launch Jupyter when you want to explore the original notebooks:

   ```bash
   jupyter notebook
   ```

## Reusable Python API

The `cgr` package exposes the core operations without requiring notebook-specific paths.

```python
from pathlib import Path

from cgr import count_kmers, fcgr_matrix, kmer_probabilities, read_fasta

sequence = read_fasta(Path("DNA_SEQUENCES/KF-430201.fasta"))
counts = count_kmers(sequence, k=5)
probabilities = kmer_probabilities(counts)
matrix = fcgr_matrix(probabilities, k=5)

print("FCGR shape:", matrix.shape)
print("Total probability:", matrix.sum())
```

For a sequence containing valid DNA bases, the FCGR matrix has shape `(2**k, 2**k)`. Ambiguous or non-DNA windows are skipped. Probabilities are normalized over the valid k-mer occurrences, so a non-empty result sums to 1.

Compare two FCGR matrices with:

```python
from cgr import fcgr_distance

distance = fcgr_distance(matrix, matrix)
print(distance)  # 0.0
```

## Run the tests

The test suite uses Python's standard `unittest` framework:

```bash
python -m unittest discover -s tests -v
```

## Notes on the research notebooks

The notebooks are retained as research artifacts and preserve the original exploratory workflow. The reusable `cgr` package is the recommended starting point for new scripts and automated tests. Notebook cells may require adapting input paths to your local checkout; use paths relative to the repository root rather than machine-specific absolute paths.

The example FASTA files and images are retained in their existing locations. Large supporting documents are also kept to avoid removing project material without a separate archival decision.

## Limitations

- This repository is an exploratory research implementation, not a clinical diagnostic tool.
- Sequence comparison results depend on the selected k-mer length and input sequences.
- The included notebooks are historical research artifacts and are not all converted into production pipelines.
- Dataset provenance and redistribution terms should be checked with the original providers before reusing the sequence data outside this project.

## Contributing

Bug reports and focused improvements are welcome. Please include the Python version, steps to reproduce, and expected versus actual behavior when reporting an issue.

## License

No license is currently specified. Unless a license is added, the repository's contents should not be assumed to be available for reuse, modification, or redistribution under an open-source license.
