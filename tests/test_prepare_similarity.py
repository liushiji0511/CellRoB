import importlib.util
from pathlib import Path
import unittest
import numpy as np

spec = importlib.util.spec_from_file_location(
    "prepare_similarity", Path(__file__).resolve().parents[1] / "scripts/prepare_similarity.py")
preparation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(preparation)


class EmbeddingOrderTests(unittest.TestCase):
    def test_reorder_and_subset(self):
        embedding = np.array([[10., 11.], [20., 21.], [30., 31.]])
        aligned = preparation.aligned_embeddings(embedding, ["a", "b", "c"], ["c", "a"])
        np.testing.assert_array_equal(aligned, [[30., 31.], [10., 11.]])

    def test_reject_ambiguous_missing_and_invalid_rows(self):
        cases = [
            (np.ones((2, 2)), ["a", "a"], ["a", "b"]),
            (np.ones((2, 2)), ["a", "b"], ["a", "a"]),
            (np.ones((2, 2)), ["a", "b"], ["a", "missing"]),
            (np.ones((3, 2)), ["a", "b"], ["a", "b"]),
            (np.array([[np.nan], [1.]]), ["a", "b"], ["a", "b"]),
        ]
        for embedding, original, requested in cases:
            with self.subTest(original=original, requested=requested):
                with self.assertRaises(ValueError):
                    preparation.aligned_embeddings(embedding, original, requested)


if __name__ == "__main__":
    unittest.main()
