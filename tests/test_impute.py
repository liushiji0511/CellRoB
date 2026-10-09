import unittest
from unittest.mock import patch
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
from cellrob import CellRoB


class ImputeTests(unittest.TestCase):
    def setUp(self):
        self.model = CellRoB.__new__(CellRoB)
        self.model.filelist = ["target", "first", "second", "cutoff", "invalid"]
        self.model.Knn = np.zeros((5, 5))
        self.model.Knn[0] = [1., 0.9, 0.7, 0.05, np.nan]
        self.model.impute_n = 15
        self.original = self.model.Knn.copy()

    def read_contact(self, path, **kwargs):
        return pd.DataFrame({"cell": [path]})

    def test_distinct_neighbors_weights_cutoff_and_no_mutation(self):
        with patch("cellrob.core.pd.read_csv", side_effect=self.read_contact) as read:
            result = self.model.impute("target")
        self.assertEqual(result["cell"].tolist(), ["target", "first", "second"])
        np.testing.assert_allclose(result["count"], [1., 0.9, 0.7])
        self.assertEqual(read.call_count, 3)
        np.testing.assert_array_equal(self.model.Knn, self.original)

    def test_neighbor_limit_and_zero_neighbors(self):
        with patch("cellrob.core.pd.read_csv", side_effect=self.read_contact):
            self.model.impute_n = 1
            self.assertEqual(self.model.impute("target")["cell"].tolist(), ["target", "first"])
            self.model.impute_n = 0
            self.assertEqual(self.model.impute("target")["cell"].tolist(), ["target"])

    def test_repeat_and_parallel_calls_are_stable(self):
        with patch("cellrob.core.pd.read_csv", side_effect=self.read_contact):
            expected = self.model.impute("target")
            with ThreadPoolExecutor(max_workers=3) as pool:
                results = list(pool.map(self.model.impute, ["target"] * 6))
        for result in results:
            pd.testing.assert_frame_equal(result, expected)
        np.testing.assert_array_equal(self.model.Knn, self.original)


if __name__ == "__main__":
    unittest.main()
