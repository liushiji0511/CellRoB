# Preparing CellRoB inputs

## Contact tables

Prepare one tab-separated file per cell. Keep one row per observed contact, with a header:

```text
chrom1	pos1	chrom2	pos2
chr6	120000	280000
chr6	160000	320000
```

Use nonnegative zero-based base-pair positions, not bin indices. If source positions are one-based, subtract one from both position columns. Match chromosome names to the reference chromosome-size file. Use the same reference assembly for positions and chromosome lengths. Include each physical contact once; CellRoB symmetrizes chromosome matrices itself. Complete mapping, quality filtering and duplicate removal before exporting contact tables. Additional columns can be retained, but CellRoB assigns a weight of one to every input row and replaces any existing `count` column.

Headerless pairs files may include read identifiers, strands and metadata. Select the chromosome and position columns according to their actual schema, rename them to the four names above, and write a header. For example, for a file whose first five fields are read ID, chromosome 1, position 1, chromosome 2, position 2:

```python
import pandas as pd
table = pd.read_csv("cell.pairs.gz", sep="\t", comment="#", header=None,
                    usecols=[1, 2, 3, 4])
table.columns = ["chrom1", "pos1", "chrom2", "pos2"]
# Subtract 1 here only if the source coordinates are one-based.
table.to_csv("cell.contacts.tsv.gz", sep="\t", index=False)
```

Plain TSV and gzip-compressed TSV are accepted. Do not expand aggregated count rows into one row without accounting for their multiplicity.

## File order and chromosome sizes

Create `filelist.txt` with one absolute contact-table path per line, using forward slashes on every platform. Use unique filenames and at least two cells. Order is preserved in every cell-by-bin output matrix.

Create a headerless, tab-separated chromosome-size file:

```text
chr6	170805979
```

The number above is an illustrative format example; use the actual length from your chosen reference assembly.

## Similarity weights from cell embeddings

The preparation workflow uses Higashi cell embeddings followed by a Scanpy neighbor graph. Obtain the cell embedding matrix from your Higashi run: shape `(number_of_cells, embedding_dimensions)`. Higashi training and its dataset-specific configuration are upstream steps, outside this repository.

Record the contact file corresponding to every embedding row in `embedding_filelist.txt`, using exactly the same path identifiers as `filelist.txt`. The preparation script reorders embedding rows to the CellRoB file order, so it can also select a subset of the embedded cells.

Install the optional preparation dependency and run:

```bash
python -m pip install scanpy
python scripts/prepare_similarity.py --embedding data/cell_embedding.npy --embedding-filelist data/embedding_filelist.txt --filelist data/filelist.txt --output data/similarity.npy --neighbors 15 --n-pcs 20 --seed 0
```

The script follows the scMicro-C preparation recipe: `sc.pp.neighbors(..., n_neighbors=15, n_pcs=20)`, followed by export of `adata.obsp['connectivities']`. It saves a dense square `.npy` matrix and a matching `.cells.txt` file. Diagonal weights are zero; CellRoB excludes the target cell during aggregation. The dense matrix requires memory proportional to the square of the number of cells.

Scanpy selects the representation according to embedding dimensionality; `n_pcs` controls principal components when PCA is used. Set `--n-pcs 0` to use the input features directly. See the [Scanpy neighbor documentation](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.neighbors.html). The earlier Dip-C preparation recipe uses the same connectivity export with Scanpy's default `n_pcs`; use the configuration matching your experiment.

`--neighbors` controls graph construction; CellRoB's `--impute-n` controls how many different neighboring cells are aggregated. These parameters serve different steps. CellRoB considers weights strictly greater than 0.05, from largest to smallest, and uses each selected weight for that cell's contacts.

For an existing similarity matrix, verify shape `(N,N)`, numeric finite nonnegative weights, and matching row/column order. Supply similarities rather than distances or integer neighbor indices. Do not substitute a different similarity construction when reproducing a specified experiment.

## Runnable example

From the repository root, after installing CellRoB:

```bash
python examples/run_example.py --output results/example
```

This generates four synthetic contact tables, a chromosome-size file, a file list and illustrative precomputed similarity weights, then executes detection, concatenation and alignment. No external dataset or Scanpy installation is required. The example uses a fictional chromosome and manually specified weights to demonstrate the input contract, rather than a biological benchmark or Higashi training result.

Outputs are written under `results/example/outputs`. The aligned matrix has four rows and 121 bins. To repeat the example, choose a new output directory.
