import tempfile
from pathlib import Path
import unittest

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


if __name__ == "__main__":
    unittest.main()
