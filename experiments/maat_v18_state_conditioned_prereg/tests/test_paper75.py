from __future__ import annotations

import csv
import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from state_conditioned_spec import (  # noqa: E402
    ResidualState,
    classify_safety,
    deterministic_split,
    hash_to_unit,
    state_supports,
)
from validate_preregistration import validate  # noqa: E402


class StateSemanticsTests(unittest.TestCase):
    def fixture(self) -> ResidualState:
        return ResidualState(
            n_vars=5,
            original_clauses=[(1, -2, 3), (-1, 2), (-3, 4, 5), (-4, -5, 2)],
        )

    def test_supports_are_bounded_and_closure_is_exact(self) -> None:
        supports = state_supports(self.fixture())
        for key in ("H", "B", "S", "V", "R_resp", "R_struct", "R_rob"):
            self.assertGreaterEqual(supports[key], 0.0)
            self.assertLessEqual(supports[key], 1.0)
        self.assertAlmostEqual(supports["R_rob"], min(supports["R_resp"], supports["R_struct"]), places=14)

    def test_assignment_and_rollback_restore_state(self) -> None:
        state = self.fixture()
        root_fingerprint = state.fingerprint()
        root_supports = state_supports(state)
        state.enqueue(1, level=1)
        state.enqueue(2, level=1, reason=(-1, 2))
        self.assertNotEqual(state.fingerprint(), root_fingerprint)
        state.backtrack(0)
        self.assertEqual(state.fingerprint(), root_fingerprint)
        self.assertEqual(state_supports(state), root_supports)

    def test_learned_clause_persists_on_rollback(self) -> None:
        state = self.fixture()
        state.enqueue(1, level=1)
        state.add_learned_clause((-1, 4))
        state.backtrack(0)
        self.assertIn((-1, 4), state.learned_clauses)
        self.assertIn((-1, 4), state.residual_clauses())

    def test_conflict_and_terminal_states_are_not_evaluated(self) -> None:
        conflict = ResidualState(1, [(1,), (-1,)])
        conflict.enqueue(1, level=0)
        with self.assertRaises(ValueError):
            state_supports(conflict)
        terminal = ResidualState(1, [(1,)])
        terminal.enqueue(1, level=0)
        with self.assertRaises(ValueError):
            state_supports(terminal)

    def test_split_key_is_deterministic(self) -> None:
        digest = "a" * 64
        self.assertEqual(deterministic_split("fixture", digest), deterministic_split("fixture", digest))
        self.assertNotEqual(deterministic_split("fixture", digest), deterministic_split("other", digest))

    def test_common_random_number_is_state_bound_not_policy_bound(self) -> None:
        state = self.fixture()
        fingerprint = state.fingerprint()
        first = hash_to_unit("fixture", 3, fingerprint)
        second = hash_to_unit("fixture", 3, fingerprint)
        null = hash_to_unit("fixture", 3, fingerprint, namespace="random_null")
        self.assertEqual(first, second)
        self.assertNotEqual(first, null)
        self.assertGreaterEqual(first, 0.0)
        self.assertLess(first, 1.0)

    def test_safety_status_distinguishes_harm_from_uncertainty(self) -> None:
        self.assertEqual(classify_safety(-0.02, 0.04), "no_harm_certified")
        self.assertEqual(classify_safety(-0.02, 0.08), "no_harm_not_established")
        self.assertEqual(classify_safety(0.06, 0.10), "harmful")


class ArtifactTests(unittest.TestCase):
    def test_protocol_and_release_gatekeeper(self) -> None:
        self.assertEqual(validate(), [])

    def test_protocol_forbids_execution_before_doi(self) -> None:
        protocol_path = ROOT / "outputs_paper75_preregistration" / "paper75_preregistration.json"
        protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
        self.assertFalse(protocol["guardrails"]["external_state_conditioned_run_before_doi"])
        self.assertEqual(protocol["status"], "preregistration_only_no_external_execution")
        self.assertIn("Zenodo DOI", protocol["release_condition"])
        self.assertEqual(protocol["evaluation_and_calibration"]["trajectory_policy"], "MOMS-only")
        self.assertEqual(protocol["frozen_constants"]["minimum_test_activation"], 0.125)
        self.assertNotIn("policy_index", protocol["policy"]["common_random_numbers"])
        self.assertIn("status_axes", protocol)

    def test_manifest_has_only_metadata_and_frozen_splits(self) -> None:
        with (ROOT / "dataset_manifest_paper75.csv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertTrue(rows)
        self.assertTrue(all(r["paper75_split"] in {"calibration", "test"} for r in rows))
        self.assertTrue(all(len(r["sha256"]) == 64 for r in rows))
        self.assertFalse(any((ROOT / r["local_path"]).exists() for r in rows))


if __name__ == "__main__":
    unittest.main()
