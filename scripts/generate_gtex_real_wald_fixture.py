#!/usr/bin/env python3
"""Extract all GTEx v8 genes for the real-data DESeq2 parity fixture."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIXTURE_DIR = ROOT / "crates/rsdeseq2/tests/data/gtex_v8_real_wald"
CACHE = ROOT / "target/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
ANNOTATIONS = ROOT / "target/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
URL = (
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
)
ANNOTATIONS_URL = (
    "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/"
    "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
)


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".partial")
    request = urllib.request.Request(url, headers={"User-Agent": "rsdeseq2-fixture-generator/1"})
    with urllib.request.urlopen(request) as response, temporary.open("wb") as output:
        shutil.copyfileobj(response, output, length=1024 * 1024)
    temporary.replace(destination)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_samples(path: Path) -> list[str]:
    with path.open(newline="") as source:
        rows = csv.DictReader(source, delimiter="\t")
        samples = [row["sample"] for row in rows]
    if len(samples) < 50 or len(set(samples)) != len(samples):
        raise ValueError(f"expected at least 50 unique sample IDs in {path}")
    return samples


def validate_sample_annotations(path: Path, samples: list[str]) -> None:
    with path.open(newline="") as source:
        annotations = {
            row["SAMPID"]: row["SMTSD"]
            for row in csv.DictReader(source, delimiter="\t")
        }
    wrong_tissue = [sample for sample in samples if annotations.get(sample) != "Kidney - Cortex"]
    if wrong_tissue:
        raise ValueError(f"samples must be exact GTEx v8 Kidney - Cortex IDs: {wrong_tissue[:5]}")


def write_design(metadata_path: Path, output_path: Path) -> None:
    with metadata_path.open(newline="") as source, output_path.open("w", newline="") as target:
        rows = csv.DictReader(source, delimiter="\t")
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "Intercept", "perm_blockmale_3", "perm_blockmale_4", "conditionB"])
        count = 0
        for row in rows:
            if row["condition"] not in {"A", "B"}:
                raise ValueError(f"unexpected condition for {row['sample']}: {row['condition']}")
            if row["perm_block"] not in {"female_4", "male_3", "male_4"}:
                raise ValueError(f"unexpected block for {row['sample']}: {row['perm_block']}")
            writer.writerow([
                row["sample"],
                1,
                int(row["perm_block"] == "male_3"),
                int(row["perm_block"] == "male_4"),
                int(row["condition"] == "B"),
            ])
            count += 1
        if count < 50:
            raise ValueError(f"expected at least 50 metadata rows, got {count}")


def extract_counts(gct_path: Path, samples: list[str], output_path: Path) -> int:
    with gzip.open(gct_path, "rt", newline="") as source:
        version = source.readline().rstrip("\r\n")
        dimensions = source.readline().rstrip("\r\n").split("\t")
        header = source.readline().rstrip("\r\n").split("\t")
        if version != "#1.2" or header[:2] != ["Name", "Description"]:
            raise ValueError("input is not the expected GTEx v8 GCT gene-read-count matrix")
        sample_indices = []
        for sample in samples:
            try:
                sample_indices.append(header.index(sample))
            except ValueError as error:
                raise ValueError(f"sample {sample} is absent from GTEx v8 GCT") from error

        rows = 0
        with output_path.open("w", newline="") as target:
            writer = csv.writer(target, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene", *samples])
            reader = csv.reader(source, delimiter="\t")
            for row in reader:
                if not row:
                    continue
                if len(row) != len(header):
                    raise ValueError(f"malformed GCT row {rows + 1}: {len(row)} fields, expected {len(header)}")
                writer.writerow([row[0], *(row[index] for index in sample_indices)])
                rows += 1

    expected_rows = int(dimensions[0])
    if rows != expected_rows:
        raise ValueError(f"GCT declares {expected_rows} genes but contains {rows}")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gct", type=Path, help="local GTEx v8 gene-reads GCT.gz; downloads it if omitted")
    parser.add_argument("--rscript", default=shutil.which("Rscript"), help="Rscript executable with DESeq2 1.52.0")
    parser.add_argument("--counts-only", action="store_true", help="extract counts without running the DESeq2 reference generator")
    args = parser.parse_args()

    gct_path = args.gct or CACHE
    if not gct_path.exists():
        print(f"Downloading GTEx v8 gene counts to {gct_path}", flush=True)
        download(URL, gct_path)
    if not ANNOTATIONS.exists():
        print(f"Downloading GTEx v8 sample annotations to {ANNOTATIONS}", flush=True)
        download(ANNOTATIONS_URL, ANNOTATIONS)

    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    metadata_path = FIXTURE_DIR / "sample_metadata.tsv"
    if not metadata_path.exists():
        raise FileNotFoundError(f"missing committed sample metadata: {metadata_path}")
    samples = read_samples(metadata_path)
    validate_sample_annotations(ANNOTATIONS, samples)
    count_path = FIXTURE_DIR / "counts.tsv"
    gene_count = extract_counts(gct_path, samples, count_path)
    write_design(metadata_path, FIXTURE_DIR / "design_matrix.tsv")

    manifest = {
        "source": "GTEx Analysis V8 gene read counts; exact Kidney - Cortex sample IDs shared with the neighboring L-BFGS-B fixture",
        "source_url": URL,
        "sample_annotations_url": ANNOTATIONS_URL,
        "source_sha256": sha256(gct_path),
        "genes": gene_count,
        "samples": len(samples),
        "tissue": "Kidney - Cortex",
        "design": "~ perm_block + condition (captured as a numeric model matrix)",
    }
    (FIXTURE_DIR / "source.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {gene_count:,} genes x {len(samples)} samples to {count_path}", flush=True)

    if args.counts_only:
        return 0
    if not args.rscript:
        raise RuntimeError("Rscript not found; activate the DESeq2 R 4.6.1 environment or pass --rscript")
    subprocess.run(
        [args.rscript, "--vanilla", str(ROOT / "scripts/generate_gtex_real_wald_reference.R"), str(FIXTURE_DIR)],
        check=True,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:  # surface a concise regeneration failure
        print(f"fixture generation failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
