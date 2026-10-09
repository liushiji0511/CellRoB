# CellRoB

CellRoB detects boundary signals in single-cell chromatin contact data and aligns nearby boundaries across cells. This repository was organized from the original `RobHiC.py` implementation; **CellRoB** is the current method name.

The workflow combines similarity-weighted contact aggregation, insulation scoring, candidate boundary scoring and cross-cell boundary alignment.

## Installation

```bash
python -m venv .venv
# Linux/macOS
source .venv/bin/activate
# Windows PowerShell: .\.venv\Scripts\Activate.ps1
python -m pip install -e .
```

Dependencies: NumPy, pandas and SciPy. See `requirements-tested.txt` for the local test environment.

## Input

- A UTF-8 text file listing one contact-table path per line. Use at least two cells for this source version. Use absolute paths with forward slashes, including on Windows. Basenames must be unique because outputs use them as cell identifiers.
- Each contact table is tab-separated with a header containing `chrom1`, `pos1`, `chrom2`, `pos2`. Each row represents one contact; an existing `count` column is overwritten. Positions are divided by the bin size using integer floor division. Supply nonnegative, zero-based positions within chromosome limits; convert one-based positions before running.
- A numerical `.npy` similarity matrix of shape `(n_cells, n_cells)`, with rows and columns in exactly the file-list order. Supply similarity weights rather than neighbor indices. See [data preparation](docs/DATA_PREPARATION.md) for the embedding-to-connectivity workflow.
- A chromosome-size table with no header: `chromosome<TAB>length_in_bp`. Chromosome names must match the contact tables.

The runner expects the zero-based coordinates described above.

## Run

For a self-contained synthetic example, run `python examples/run_example.py --output results/example`. See [data preparation and examples](docs/DATA_PREPARATION.md) to prepare your own contacts and similarity matrix.

```bash
cellrob --filelist data/filelist.txt --knn data/similarity.npy \
  --chrom-sizes data/genome.chrom.sizes --chroms chr6 \
  --output results/cellrob --binsize 40000 --window 15 \
  --impute-n 15 --workers 1 --seed 0
```

In PowerShell, enter this command on one line. `python -m cellrob` is equivalent to `cellrob`.

Stages run in order: `run`, `concatenate`, `align`. Use `--stage` to run a single stage with its prerequisite outputs already present. Use a fresh output directory for each configuration. Use `--workers 1` and a fixed `--seed` for reproducible random perturbations.

The example uses a window of 15 bins at 40 kb. The core default is 5 bins. Set the window for your experimental configuration. Read [the method description](docs/METHOD.md) for the scoring and alignment procedure.

## Output

| Path under the output directory | Meaning |
| --- | --- |
| `<chrom>/<input_basename>.boundary.npy` | Boundary score per genomic bin; not a binary call |
| `<chrom>/<input_basename>.insulation.npy` | Insulation score per genomic bin |
| `concatenate/<chrom>.boundary.npy` | Cell-by-bin boundary scores, in file-list order |
| `concatenate/<chrom>.insulation.npy` | Cell-by-bin insulation scores |
| `align/<chrom>.boundary.npy` | Aligned cell-by-bin boundary strengths; zeros indicate absence |

Array length is `floor(chromosome_length / binsize) + 1`, as in the source. Bin `j` corresponds to nominal coordinates `[j*binsize, (j+1)*binsize)`; the final bin may be empty or extend beyond the chromosome. Alignment retains each participating cell's maximum insulation score in its region at a shared anchor.

## Repository

```text
src/cellrob/       Core implementation and command-line runner
docs/             Method explanation
examples/         Runnable synthetic dataset example
scripts/          Similarity matrix preparation
tests/            Neighbor aggregation regression tests
pyproject.toml    Package metadata
requirements.txt Runtime dependencies
```

The local `_source_snapshot/` contains byte-preserved candidate sources and SHA-256 hashes for comparison. It is excluded from Git: historical scripts include machine-specific paths and unrelated experimental code. Dataset processing, benchmark methods, figures, result arrays and editor files are excluded from the publication package.

Run the neighbor aggregation tests with `python -m unittest discover -s tests -v`.
