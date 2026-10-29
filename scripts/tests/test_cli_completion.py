"""A completion record must reject partially published or changed generations."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from check_cli_completion import verify


class CompletionChecks(unittest.TestCase):
    def test_valid_generation_then_mutated_missing_and_duplicate_output(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            output = root/'results.tsv'
            output.write_bytes(b'gene\tpvalue\ng1\t0.01\n')
            item = dict(path=str(output.resolve()), bytes=output.stat().st_size,
                        sha256=hashlib.sha256(output.read_bytes()).hexdigest())
            manifest = dict(schema_version=1, generation='checked-generation', files=[item])
            marker = root/'results.tsv.complete.json'
            marker.write_text(json.dumps(manifest))
            self.assertEqual(verify(marker), manifest)
            output.write_bytes(b'partial')
            with self.assertRaisesRegex(ValueError, 'changed output'):
                verify(marker)
            output.unlink()
            with self.assertRaises(FileNotFoundError):
                verify(marker)
            output.write_bytes(b'gene\tpvalue\ng1\t0.01\n')
            manifest['files'].append(item)
            marker.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, 'unique'):
                verify(marker)


if __name__ == '__main__':
    unittest.main()
