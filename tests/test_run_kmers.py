import tempfile
from pathlib import Path
import unittest
from unittest import mock
import subprocess

from mbcclr_utils import runners_utils


class RunKmersTests(unittest.TestCase):
    def test_run_kmers_writes_profiles_with_expected_shape(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reads = tmp / "reads.fasta"
            output = tmp / "output"
            output.mkdir()

            reads.write_text(
                ">r1\nACGTACGT\n>r2\nACGTNNNN\n>r3\nAA\n",
                encoding="utf-8",
            )

            runners_utils.run_kmers(str(reads), str(output), 3, 1)

            lines = (output / "profiles" / "3mers").read_text(encoding="utf-8").strip().splitlines()
            self.assertEqual(len(lines), 3)

            v1 = [float(x) for x in lines[0].split()]
            v2 = [float(x) for x in lines[1].split()]
            v3 = [float(x) for x in lines[2].split()]

            self.assertEqual(len(v1), 32)
            self.assertEqual(len(v2), 32)
            self.assertEqual(len(v3), 32)
            self.assertAlmostEqual(sum(v1), 1.0, places=4)
            self.assertAlmostEqual(sum(v2), 1.0, places=4)
            self.assertAlmostEqual(sum(v3), 0.0, places=4)

    def test_run_15mer_vecs_uses_cov_and_writes_expected_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            reads = tmp / "reads.fasta"
            output = tmp / "output"
            output.mkdir()
            reads.write_text(">r1\nACGTACGTACGTACGT\n", encoding="utf-8")

            seen_cmd = {}

            def fake_run(cmd):
                seen_cmd["cmd"] = cmd
                out_dir = cmd[cmd.index("-o") + 1]
                Path(out_dir, "kmers.counts").write_text("1\t2\n", encoding="utf-8")
                Path(out_dir, "kmers.vectors").write_text("0.1 0.9\n", encoding="utf-8")
                return subprocess.CompletedProcess(cmd, 0)

            with mock.patch("mbcclr_utils.runners_utils.shutil.which", return_value="/usr/bin/kmertools"):
                with mock.patch("mbcclr_utils.runners_utils.subprocess.run", side_effect=fake_run):
                    runners_utils.run_15mer_vecs(str(reads), str(output), 10, 32, 1)

            self.assertIn("cmd", seen_cmd)
            self.assertEqual(seen_cmd["cmd"][1], "cov")
            self.assertEqual((output / "profiles" / "15mers-counts").read_text(encoding="utf-8"), "1\t2\n")
            self.assertEqual((output / "profiles" / "15mers").read_text(encoding="utf-8"), "0.1 0.9\n")


if __name__ == "__main__":
    unittest.main()
