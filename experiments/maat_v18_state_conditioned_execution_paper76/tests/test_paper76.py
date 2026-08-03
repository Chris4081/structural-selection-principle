from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from paper76_core import (  # noqa: E402
    IncrementalResidualTracker,
    hash_to_unit,
    max_snapshot_difference,
    weighted_quantile_lower,
)
from paper76_execution import (  # noqa: E402
    PAPER75_DOI,
    PAPER75_MANIFEST_SHA256,
    PAPER75_PROTOCOL_SHA256,
    Instance,
    run_solver,
    validate_release,
)


class IncrementalTrackerTests(unittest.TestCase):
    def fixture(self) -> IncrementalResidualTracker:
        return IncrementalResidualTracker([(1, -2, 3), (-1, 2), (2, 3, -4), (-3, 4)], 4)

    def assert_audit(self, tracker: IncrementalResidualTracker) -> None:
        vsids = {v: 1.0 + 0.1 * v for v in range(1, 5)}
        self.assertLessEqual(max_snapshot_difference(tracker.snapshot(vsids), tracker.full_snapshot(vsids)), 1.0e-14)

    def test_enqueue_and_rollback_match_full_recomputation(self) -> None:
        tracker = self.fixture()
        self.assert_audit(tracker)
        tracker.assign(1, True)
        self.assert_audit(tracker)
        tracker.assign(2, True)
        self.assert_audit(tracker)
        tracker.unassign(2)
        tracker.unassign(1)
        self.assert_audit(tracker)

    def test_learned_clause_persists(self) -> None:
        tracker = self.fixture()
        tracker.assign(1, True)
        tracker.add_clause((-1, 4))
        tracker.unassign(1)
        self.assertIn((-1, 4), tracker.clauses)
        self.assert_audit(tracker)

    def test_common_random_number_has_no_policy_identity(self) -> None:
        first = hash_to_unit("fixture", 7, "abc")
        second = hash_to_unit("fixture", 7, "abc")
        null = hash_to_unit("fixture", 7, "abc", namespace="random_null")
        self.assertEqual(first, second)
        self.assertNotEqual(first, null)

    def test_lower_weighted_quantile(self) -> None:
        self.assertEqual(weighted_quantile_lower([0.1, 0.2, 0.3, 0.4], [0.25] * 4, 0.75), 0.3)


class ReleaseAndSyntheticSmokeTests(unittest.TestCase):
    def test_paper75_release_and_hashes_are_frozen(self) -> None:
        self.assertEqual(validate_release(), [])
        release = json.loads((ROOT / "paper75_release_verification.json").read_text(encoding="utf-8"))
        self.assertEqual(release["doi"], PAPER75_DOI)
        self.assertEqual(release["protocol_sha256"], PAPER75_PROTOCOL_SHA256)
        self.assertEqual(release["manifest_sha256"], PAPER75_MANIFEST_SHA256)

    def test_synthetic_moms_trace_only(self) -> None:
        instance = Instance(
            dataset_id="synthetic_excluded_fixture",
            family="synthetic_excluded",
            local_path="not_external.cnf",
            sha256="0" * 64,
            split="calibration",
            instance_weight=1.0,
            n_vars=4,
            n_clauses=5,
            clauses=[(1, 2), (-1, 3), (-2, -3), (2, 4), (-4,)],
        )
        record, trace = run_solver(instance, "moms", None, collect_trace=True)
        self.assertEqual(record.policy, "moms")
        self.assertLessEqual(record.audit_max_abs_error, 1.0e-12)
        self.assertGreaterEqual(len(trace), 0)


if __name__ == "__main__":
    unittest.main()

