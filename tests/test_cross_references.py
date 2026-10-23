"""Cross-reference guards for the v11 article and its supplement.

* Numbered results (lemmas, definitions, propositions, corollaries, theorems) are cited
  through labels, never by a typed number, in the article, the supplement and the tables.
* Every reference of the supplement to the main article (prefix ``M-``, read with xr-hyper)
  names a label that exists in the main article.
* The main article cannot read the supplement's labels (it must compile on its own), so its
  typed references "Supplementary Material S<n>" are checked here against the order of the
  supplement's sections.
"""
import re
import unittest
from pathlib import Path

MS = Path(__file__).resolve().parents[1] / "manuscript"
MAIN = (MS / "manuscript.tex").read_text(encoding="utf-8")
SUPP = (MS / "supplementary.tex").read_text(encoding="utf-8")
TABLES = {p.name: p.read_text(encoding="utf-8") for p in MS.glob("table_*.tex")}

TYPED_RESULT = re.compile(r"\b(?:Theorem|Proposition|Corollary|Lemma|Definition)s?(?:~|\s)[0-9]")


def labels(text):
    return (set(re.findall(r"\\label\{([^}]+)\}", text))
            | set(re.findall(r"\\stmt\{[^}]+\}\{([^}]+)\}", text)))


class CrossReferences(unittest.TestCase):
    def test_numbered_results_are_cited_by_label(self):
        main_tables = "".join(TABLES.values())
        self.assertEqual(TYPED_RESULT.findall(MAIN), [])
        self.assertEqual(TYPED_RESULT.findall(main_tables), [])
        # research questions are defined with \rqdef and cited with \rqref
        self.assertEqual(re.findall(r"\bRQ[0-9]", MAIN), [])
        self.assertEqual(len(re.findall(r"\\item\[\\textbf\{\\rqdef\{", MAIN)), 3)
        # the supplement cites one theorem of another paper by its number
        self.assertEqual(TYPED_RESULT.findall(SUPP.replace("whose Theorem~4", "")), [])

    def test_supplement_references_to_the_main_article_resolve(self):
        main_labels = labels(MAIN) | labels("".join(TABLES.values()))
        used = set(re.findall(r"\\ref\*?\{M-([^}]+)\}", SUPP + "".join(TABLES.values())))
        self.assertGreater(len(used), 15)
        self.assertEqual(used - main_labels, set())
        self.assertNotRegex(SUPP, r"(?:Section|Table|Figure)s?~[0-9][0-9.]*(?:~and~[0-9.]+)? of the main article")

    def test_typed_supplement_section_numbers_in_the_main_article(self):
        order = re.findall(r"\\section\{[^}]*\}\s*\\label\{(sec:supp_[a-z_]+)\}", SUPP)
        number = {lab: i + 1 for i, lab in enumerate(order)}
        self.assertEqual(len(order), 7)
        expected = {  # phrase in the main article -> supplement section it must point to
            "given in Supplementary Material S1": "sec:supp_trace",
            "documented in Supplementary Material S2": "sec:supp_pipeline",
            "listed in Supplementary Material S6": "sec:supp_rts",
            "Supplementary Material S6 plots the minimum cost": "sec:supp_rts",
            "Supplementary Material Section~S4 rests on": "sec:supp_pharma",
            "Supplementary Material S6 maps each reading": "sec:supp_rts",
            "Supplementary Material S3 records the verification": "sec:supp_oracle",
            "Supplementary Material S7 gives these results": "sec:supp_baselines",
        }
        flat = re.sub(r"\s+", " ", MAIN)
        for phrase, lab in expected.items():
            self.assertIn(phrase, flat)
            self.assertIn(f"S{number[lab]}", phrase, phrase)
        # the pair "S4 and S5" always names the two literature-derived instances
        for m in re.finditer(r"Supplementary Material S(\d) and S(\d)", flat):
            self.assertEqual((int(m.group(1)), int(m.group(2))),
                             (number["sec:supp_pharma"], number["sec:supp_fyffe"]))
        cited = set()
        for n in re.findall(r"Supplementary Material (?:Section~)?S(\d)", flat):
            self.assertLessEqual(int(n), len(order))
            cited.add(int(n))
        for a, b in re.findall(r"Supplementary Material S(\d) and S(\d)", flat):
            cited |= {int(a), int(b)}
        # the supplement introduction says its sections are cited from the main text as S1--S7
        self.assertEqual(cited, set(range(1, len(order) + 1)))

    def test_every_float_is_cited(self):
        for name, text in (("main", MAIN + "".join(TABLES.values())), ("supplement", SUPP + "".join(TABLES.values()))):
            src = MAIN if name == "main" else SUPP
            for env in ("table", "figure", "algorithm"):
                for m in re.finditer(r"\\begin\{%s\}(.*?)\\end\{%s\}" % (env, env), src, re.S):
                    labs = re.findall(r"\\label\{([^}]+)\}", m.group(1))
                    self.assertTrue(labs, f"unlabelled {env} in {name}")
                    rest = src[:m.start()] + src[m.end():]
                    self.assertRegex(rest, r"\\ref\*?\{%s\}" % re.escape(labs[0]), f"{labs[0]} not cited in {name}")
            for p in [p for p in TABLES if ("\\input{%s}" % p) in src]:
                lab = re.search(r"\\label\{([^}]+)\}", TABLES[p]).group(1)
                self.assertRegex(src, r"\\ref\*?\{%s\}" % re.escape(lab), f"{lab} not cited in {name}")


if __name__ == "__main__":
    unittest.main()
