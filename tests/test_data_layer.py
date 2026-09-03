import unittest

import numpy as np
import pandas as pd

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from data_layer import add_canonical_group, compute_plasma_status, earliest_row


class DataLayerTests(unittest.TestCase):
    def test_earliest_row_keeps_one_source_record(self):
        frame = pd.DataFrame({
            "RID": [1, 1],
            "DATE": ["2020-01-01", "2020-02-01"],
            "A": [np.nan, 99],
            "B": ["early", "late"],
        })
        result = earliest_row(frame, "DATE", ["A", "B"])
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["B"], "early")

    def test_missing_marker_is_not_negative_by_default(self):
        frame = pd.DataFrame({"PTAU217_bl": [np.nan, 0.2, 0.1]})
        status = compute_plasma_status(frame)
        self.assertTrue(pd.isna(status.iloc[0]))
        self.assertEqual(int(status.iloc[1]), 1)
        self.assertEqual(int(status.iloc[2]), 0)

    def test_canonical_group_codes(self):
        frame = pd.DataFrame({"PET_STATUS": [0, 0, 1, 1], "PTAU217_bl": [0.1, 0.2, 0.1, 0.2]})
        result = add_canonical_group(frame)
        self.assertEqual(result["GROUP"].tolist(), [0, 1, 2, 3])


if __name__ == "__main__":
    unittest.main()
