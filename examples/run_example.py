"""Generate synthetic contacts and run all CellRoB stages without external data."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from cellrob import CellRoB


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="results/example")
    args = parser.parse_args()
    root = Path(args.output).resolve()
    if root.exists():
        parser.error("Choose a new output directory to generate a fresh example")
    inputs = root / "inputs"
    inputs.mkdir(parents=True)
    rng = np.random.default_rng(7)
    files = []
    for cell in range(4):
        rows = []
        for i in range(120):
            for j in range(i + 1, min(i + 16, 120)):
                same_block = i // 20 == j // 20
                probability = 0.8 if same_block else 0.25
                if j == i + 1 or rng.random() < probability:
                    rows.append(("chrDemo", i * 40000, "chrDemo", j * 40000))
        path = inputs / f"cell_{cell + 1}.contacts.tsv"
        pd.DataFrame(rows, columns=["chrom1", "pos1", "chrom2", "pos2"]).to_csv(path, sep="\t", index=False)
        files.append(path.as_posix())
    filelist = inputs / "filelist.txt"
    filelist.write_text("\n".join(files) + "\n", encoding="utf-8")
    chromsizes = inputs / "demo.chrom.sizes"
    chromsizes.write_text("chrDemo\t4800000\n", encoding="utf-8")
    # Illustrative weights for the tutorial; biological weights use prepare_similarity.py.
    weights = np.array([[0, .9, .7, .4], [.9, 0, .6, .5], [.7, .6, 0, .8], [.4, .5, .8, 0]])
    similarity = inputs / "similarity.npy"
    np.save(similarity, weights)
    np.random.seed(0)
    model = CellRoB(filelist_dir=str(filelist), Knn_path=str(similarity),
                    chromsizes_dir=str(chromsizes), chroms=["chrDemo"],
                    result_dir=str(root / "outputs"), binsize=40000,
                    window=5, impute_n=2, worker_n=1)
    model.run()
    model.concatenate_boundary_insulationScore()
    model.align()
    arrays = list((root / "outputs").rglob("*.npy"))
    if len(arrays) != 11:
        raise RuntimeError(f"Expected 11 output arrays; found {len(arrays)}")
    aligned = np.load(root / "outputs/align/chrDemo.boundary.npy")
    assert aligned.shape == (4, 121)
    print(f"Example completed: {root}; aligned matrix shape: {aligned.shape}")


if __name__ == "__main__":
    main()
