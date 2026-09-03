import unittest

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from environment import ExplorationEnvironment


class EnvironmentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.env = ExplorationEnvironment(seed=2026)

    def test_all_actions_have_structured_feedback(self):
        actions = [
            ("define_discordance", {"missing_policy": "exclude"}),
            ("discover_subtypes", {"method": "kmeans", "k": 2}),
            ("select_slice", {"diagnosis": "CN"}),
            ("test_confounder", {"control_vars": ["CARDIO", "ENDO"]}),
            ("sensitivity_analysis", {"outcome": "CDRSB", "window": ">=2yr_followup"}),
            ("profile_mechanism", {"target_group": "PET−/Plasma+", "alignment_window": "180d"}),
        ]
        for action, params in actions:
            feedback = self.env.act(action, params)
            self.assertNotIn("error", feedback, action)
            self.assertIn("metrics", feedback, action)

    def test_invalid_mechanism_window_is_recorded_as_error(self):
        feedback = self.env.act("profile_mechanism", {"alignment_window": "7d"})
        self.assertIn("error", feedback)


if __name__ == "__main__":
    unittest.main()
