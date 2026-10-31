"""The optimization trace of Supplementary Material S1 matches its generating script."""
import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUPP = (ROOT / "manuscript" / "supplementary.tex").read_text(encoding="utf-8")


def _load():
    spec = importlib.util.spec_from_file_location("supp_trace_v11", ROOT / "experiments" / "supp_trace_v11.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.main()


class SuppTraceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rep = _load()
        a = SUPP.index(r"\subsection{Optimization Trace")
        cls.text = SUPP[a:SUPP.index(r"\section{Computational Pipeline")]

    def test_result(self):
        self.assertEqual(self.rep["seed_cost"], 40)
        self.assertEqual(self.rep["Z_star"], 40)
        self.assertEqual(self.rep["x_star"], [1, 1, 1])
        self.assertAlmostEqual(self.rep["R_star"], 0.974, places=12)
        self.assertAlmostEqual(self.rep["R_designs"]["[2,1,1]"], 0.9784, places=12)
        self.assertIn(r"([1,1,1],\,40,\,0.974)", self.text)

    def test_phase_one_events(self):
        p1 = [e for e in self.rep["log"] if e.get("phase") == 1 and "node" in e]
        first, pruned = p1[0], p1[-1]
        self.assertEqual((first["node"], first["Z_next"], first["Z_suf"], first["event"]), ([1], 10, 30, "branch"))
        self.assertAlmostEqual(first["q_next"], 4.44)
        self.assertEqual((pruned["node"], pruned["Z_next"], pruned["event"]), ([2], 20, "pruned by Theorem 2"))
        self.assertAlmostEqual(pruned["q_next"], 3.775)
        for s in ("q_{\\mathrm{next}}=9.0-4.56=4.44", "\\Zsuf_1(4.44)=40", "q_{\\mathrm{next}}=9.0-5.225=3.775",
                  "20+\\Zsuf_1(3.775)=50>40", "\\Rsys^{0.5,0}([2,1,1])=0.9784", "W^{\\mathrm{max}}_1=9.5"):
            self.assertIn(s, self.text)

    def test_numbers_consistent_with_effective_values(self):
        self.assertEqual(self.rep["k_alpha_gamma"], 9.0)
        self.assertEqual(self.rep["effective_capacities"], [[4.56, 5.225], [4.75], [4.75]])
        self.assertEqual(self.rep["W_max_suffix"][0], 9.5)


if __name__ == "__main__":
    unittest.main()
