"""scripts/check_report_numbers.py: the real reports pass, and planted mismatches fail."""

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("crn", ROOT / "scripts" / "check_report_numbers.py")
crn = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(crn)


def _write(tmp_path, md, tex, allow=""):
    (tmp_path / "a.md").write_text(md)
    (tmp_path / "a.tex").write_text(tex)
    (tmp_path / "allow.txt").write_text(allow)
    return tmp_path / "a.md", [tmp_path / "a.tex"], tmp_path / "allow.txt"


def test_repository_reports_pass():
    assert crn.main() == 0


def test_matching_numbers_pass(tmp_path):
    assert crn.check(*_write(tmp_path, "depth 0.30 mK at 1.4 x 10^-8", r"depth \SI{0.3}{mK} at $1.4\times10^{-8}$")) == 0


def test_planted_mismatch_fails(tmp_path):
    # Negative control: a transcription error (0.31 for 0.30) must be caught.
    assert crn.check(*_write(tmp_path, "depth 0.30 mK", r"depth \SI{0.31}{mK}")) == 1


def test_dropped_minus_sign_fails(tmp_path):
    assert crn.check(*_write(tmp_path, "shift -0.18 um", r"shift \SI{0.18}{um}")) == 1


def test_allowlisted_number_passes_and_stale_entry_fails(tmp_path):
    allow = "a.tex | 25.4 | 1 inch in the Markdown\n"
    assert crn.check(*_write(tmp_path, "a 1 inch mirror", r"a \SI{25.4}{mm} mirror", allow)) == 0
    assert crn.check(*_write(tmp_path, "a 1 inch mirror", r"a one-inch mirror", allow)) == 1
