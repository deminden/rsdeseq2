#!/usr/bin/env python3
"""Record the identity of an archive after its target's release tests pass."""
import argparse
import json
from pathlib import Path
import subprocess
from release_artifact_manifest import sha256, source_identity
from check_release_target import RELEASE_TESTS

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--binary', type=Path, required=True)
parser.add_argument('--archive', type=Path, required=True)
parser.add_argument('--target', required=True)
parser.add_argument('--validation', type=Path, default=Path('validated/validation.json'))
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
validation = json.loads(args.validation.read_text())
if validation['source_commit'] != source_identity():
    raise SystemExit('Binary and checked source revisions differ')
args.output.write_text(json.dumps(dict(source_commit=source_identity(), target=args.target,
    validation_sha256=sha256(args.validation), binary_sha256=sha256(args.binary),
    archive_sha256=sha256(args.archive), rustc=subprocess.check_output(['rustc','-vV'],text=True),
    profile='release', numerical_tests=list(RELEASE_TESTS)), indent=2)+'\n')
