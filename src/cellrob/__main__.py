import argparse
import random
import numpy as np
from . import CellRoB


def main():
    parser = argparse.ArgumentParser(description="CellRoB boundary detection and alignment")
    parser.add_argument("--filelist", required=True)
    parser.add_argument("--knn", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--chroms", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--binsize", type=int, default=40000)
    parser.add_argument("--window", type=int, default=5)
    parser.add_argument("--impute-n", type=int, default=15)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--count-threshold", type=int, default=0)
    parser.add_argument("--stage", choices=["all", "run", "concatenate", "align"], default="all")
    args = parser.parse_args()
    if args.binsize <= 0 or args.binsize > 40000:
        parser.error("binsize must be between 1 and 40000 for this source implementation")
    if args.window < 2 or args.workers < 1 or args.impute_n < 0:
        parser.error("window >= 2, workers >= 1 and impute-n >= 0 are required")
    if args.count_threshold < 0:
        parser.error("count-threshold must be nonnegative")
    random.seed(args.seed)
    np.random.seed(args.seed)
    model = CellRoB(filelist_dir=args.filelist, Knn_path=args.knn,
                    chromsizes_dir=args.chrom_sizes, chroms=args.chroms,
                    result_dir=args.output, binsize=args.binsize,
                    window=args.window, impute_n=args.impute_n, worker_n=args.workers)
    if args.stage in ("all", "run"):
        model.run()
    if args.stage in ("all", "concatenate"):
        model.concatenate_boundary_insulationScore()
    if args.stage in ("all", "align"):
        model.align(count_threshold=args.count_threshold)


if __name__ == "__main__":
    main()
