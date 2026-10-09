"""Build CellRoB weights from cell embeddings using the original Scanpy recipe."""
import argparse
from pathlib import Path
import numpy as np


def aligned_embeddings(embedding, embedding_files, contact_files):
    if embedding.ndim != 2 or embedding.shape[0] != len(embedding_files):
        raise ValueError("Embedding rows must match embedding-filelist entries")
    if not np.isfinite(embedding).all():
        raise ValueError("Embedding must contain finite numeric values")
    if len(set(embedding_files)) != len(embedding_files):
        raise ValueError("Embedding file identifiers must be unique")
    if len(set(contact_files)) != len(contact_files):
        raise ValueError("Contact file identifiers must be unique")
    lookup = {name: i for i, name in enumerate(embedding_files)}
    missing = set(contact_files) - lookup.keys()
    if missing:
        raise ValueError(f"Contact files absent from embedding-filelist: {sorted(missing)}")
    return embedding[[lookup[name] for name in contact_files]]


def read_files(path):
    return [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embedding", required=True)
    parser.add_argument("--embedding-filelist", required=True)
    parser.add_argument("--filelist", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--neighbors", type=int, default=15)
    parser.add_argument("--n-pcs", type=int, default=20)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    files = read_files(args.filelist)
    embedding = aligned_embeddings(np.load(args.embedding, allow_pickle=False),
                                   read_files(args.embedding_filelist), files)
    if not 2 <= args.neighbors < len(files):
        parser.error("neighbors must be at least 2 and smaller than the number of cells")
    if args.n_pcs < 0:
        parser.error("n-pcs must be nonnegative")
    import scanpy as sc
    import inspect
    adata = sc.AnnData(embedding)
    kwargs = dict(n_neighbors=args.neighbors, n_pcs=args.n_pcs)
    seed_parameter = "rng" if "rng" in inspect.signature(sc.pp.neighbors).parameters else "random_state"
    kwargs[seed_parameter] = args.seed
    sc.pp.neighbors(adata, **kwargs)
    weights = adata.obsp["connectivities"].toarray()
    np.fill_diagonal(weights, 0)
    output = Path(args.output)
    if output.suffix != ".npy":
        parser.error("output must end in .npy")
    output.parent.mkdir(parents=True, exist_ok=True)
    np.save(output, weights)
    output.with_suffix(".cells.txt").write_text("\n".join(files) + "\n", encoding="utf-8")
    print(f"Saved {weights.shape} weights and ordered cell identifiers to {output}")


if __name__ == "__main__":
    main()
