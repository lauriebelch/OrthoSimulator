#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import gzip
import re
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

from Bio import SeqIO


def parse_args():
    script_dir = Path(__file__).resolve().parent
    # scripts/PFAM contains code only. Defaults point to the top-level PFAM
    # data/results folder when this script is run manually from the normal
    # OrthoSimulator layout. scripts/PFAM/pfam_wrapper.py still passes these explicitly.
    project_root = script_dir.parents[1]
    default_pfam_dir = project_root / "PFAM"
    default_scan = default_pfam_dir / "pfamscan_results"
    default_results = default_pfam_dir / "pfam_results"

    parser = argparse.ArgumentParser(
        description="Extract PFAM substitution models from pfam_scan.pl outputs."
    )
    parser.add_argument(
        "--pfam-scan-folder",
        default=str(default_scan),
        help="Folder containing pfam_scan.pl .txt output files.",
    )
    parser.add_argument(
        "--results-dir",
        default=str(default_results),
        help="PFAM model results directory.",
    )
    parser.add_argument(
        "--threads",
        default="1",
        help="Threads to pass to IQ-TREE.",
    )
    return parser.parse_args()


def safe_int_threads(value):
    try:
        threads = int(value)
    except ValueError:
        print(f"ERROR: --threads must be an integer, got: {value}")
        sys.exit(1)
    if threads < 1:
        print(f"ERROR: --threads must be >= 1, got: {value}")
        sys.exit(1)
    return threads


def iqtree_executable():
    """IQ-TREE 3 as shipped in OrthoSimulator's bin/ folder (the same one the
    simulation uses), or iqtree3 on the PATH if bin/ has none."""
    bundled = Path(__file__).resolve().parents[2] / "bin" / "iqtree3"
    if bundled.is_file():
        return str(bundled)
    found = shutil.which("iqtree3")
    if found:
        return found
    raise FileNotFoundError("IQ-TREE 3 not found: expected bin/iqtree3 in the OrthoSimulator folder or iqtree3 on the PATH")


class PfamModelExtractor:
    def __init__(self, pfam_scan_folder, results_dir, threads):
        self.pfam_scan_folder = Path(pfam_scan_folder).expanduser().resolve()
        self.results_dir = Path(results_dir).expanduser().resolve()
        self.threads = threads

        self.alignment_dir = self.results_dir / "alignments"
        self.trimmed_dir = self.results_dir / "trimmed_alignments"
        self.iqtree_dir = self.results_dir / "iqtree_files"
        self.csv_dir = self.results_dir / "csv_files"
        self.species_model_dir = self.csv_dir / "species_models"

    def make_dirs(self):
        folders = [
            self.alignment_dir,
            self.trimmed_dir,
            self.iqtree_dir,
            self.csv_dir,
            self.species_model_dir,
        ]
        for folder in folders:
            folder.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def species_name_from_pfam_file(path):
        name = path.name
        if name.endswith("_pfam.txt"):
            return name.replace("_pfam.txt", "")
        return path.stem

    @staticmethod
    def parse_pfam_scan_file(path):
        """
        Parse one pfam_scan.pl output file.

        Keeps rows where column 8, 'type', is 'Domain'.  Pulls PFAM accessions
        from column 6, 'hmm acc', and strips accession versions:
            PF01549.10 -> PF01549
        """
        pfams = set()
        with open(path) as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 8:
                    continue
                hmm_acc = parts[5]
                hit_type = parts[7]
                if hit_type != "Domain":
                    continue
                pfam_id = hmm_acc.split(".")[0]
                pfams.add(pfam_id)
        return sorted(pfams)

    def collect_species_pfams(self):
        species_pfams = {}
        all_pfams = set()
        files = sorted(self.pfam_scan_folder.glob("*.txt"))
        if not files:
            print(f"WARNING: no .txt files found in {self.pfam_scan_folder}")
        for path in files:
            species = self.species_name_from_pfam_file(path)
            print(f"Reading {path.name}")
            pfams = self.parse_pfam_scan_file(path)
            species_pfams[species] = pfams
            all_pfams.update(pfams)
        return species_pfams, sorted(all_pfams)

    def download_seed_alignment(self, pfam_id):
        gz_file = self.alignment_dir / f"{pfam_id}.seed.gz"
        sth_file = self.alignment_dir / f"{pfam_id}.seed.sth"
        fasta_file = self.alignment_dir / f"{pfam_id}.seed.fasta"
        if fasta_file.exists():
            return fasta_file
        if not sth_file.exists():
            url = (
                f"https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/"
                f"{pfam_id}/?annotation=alignment:seed&download"
            )
            print(f"Downloading {pfam_id}")
            try:
                urllib.request.urlretrieve(url, gz_file)
                with gzip.open(gz_file, "rb") as src:
                    with open(sth_file, "wb") as dst:
                        shutil.copyfileobj(src, dst)
            except Exception as err:
                print(f"WARNING: failed to download {pfam_id}: {err}")
                return None
        print(f"Converting {pfam_id} Stockholm to FASTA")
        try:
            records = SeqIO.parse(str(sth_file), "stockholm")
            count = SeqIO.write(records, str(fasta_file), "fasta")
        except Exception as err:
            print(f"WARNING: failed to convert {pfam_id} to FASTA: {err}")
            return None
        if count == 0:
            print(f"WARNING: no FASTA records written for {pfam_id}")
            return None
        return fasta_file

    def run_trimal(self, fasta_file):
        pfam_id = fasta_file.name.split(".")[0]
        trimmed_file = self.trimmed_dir / f"{pfam_id}.trimmed.fasta"
        if trimmed_file.exists():
            return trimmed_file
        cmd = ["trimal", "-in", str(fasta_file), "-out", str(trimmed_file), "-gappyout"]
        print(f"Trimming {pfam_id}")
        subprocess.run(cmd, check=True)
        return trimmed_file

    def run_iqtree(self, alignment_file):
        pfam_id = alignment_file.name.split(".")[0]
        output_prefix = self.iqtree_dir / pfam_id
        iqtree_file = Path(str(output_prefix) + ".iqtree")
        if iqtree_file.exists():
            return iqtree_file
        cmd = [
            iqtree_executable(),
            "--quiet",
            "-s", str(alignment_file),
            "-mset", "LG+I,JTT+G",
            "-T", str(self.threads),
            "--cmax", "4",
            "--fast",
            "--prefix", str(output_prefix),
        ]
        print(f"Running IQ-TREE for {pfam_id}")
        subprocess.run(cmd, check=True)
        return iqtree_file

    def process_one_pfam(self, pfam_id):
        fasta_file = self.download_seed_alignment(pfam_id)
        if fasta_file is None:
            return None
        # Keep current behaviour: IQ-TREE is run on the seed FASTA directly.
        # TrimAl support remains available via run_trimal() if we choose to enable it later.
        iqtree_file = self.run_iqtree(fasta_file)
        return iqtree_file

    @staticmethod
    def parse_iqtree(path):
        with open(path) as handle:
            text = handle.read()
        model = re.search(r"Model of substitution:\s*(\S+)", text)
        pinv = re.search(r"Proportion of invariable sites:\s*([\d.]+)", text)
        alpha = re.search(r"Gamma shape alpha:\s*([\d.]+)", text)
        rates = re.search(r"Site proportion and rates:\s*(.+)", text)
        model_str = model.group(1) if model else ""
        I_str = pinv.group(1) if pinv else ""
        alpha_str = alpha.group(1) if alpha else ""
        rates_str = rates.group(1) if rates else ""
        toks = model_str.split("+")
        keep = [t for t in toks if not re.match(r"^(I|G\d+|R\d+)$", t)]
        base = "+".join(keep)
        parts = [base]
        if I_str:
            parts.append(f"I{{{I_str}}}")
        if alpha_str:
            gamma_match = re.search(r"G(\d+)", model_str)
            if gamma_match:
                k = gamma_match.group(1)
                parts.append(f"G{k}{{{alpha_str}}}")
        if rates_str:
            pairs = re.findall(r"\(([\d.]+),([\d.]+)\)", rates_str)
            if pairs:
                k = len(pairs)
                flat = ",".join([",".join(pair) for pair in pairs])
                parts.append(f"R{k}" + "{" + flat + "}")
        alisim_model = "+".join([p for p in parts if p])
        pfam_id = path.name.replace(".iqtree", "")
        return {
            "pfam_id": pfam_id,
            "model": model_str,
            "I": I_str,
            "alpha": alpha_str,
            "rates": rates_str,
            "alisim_model": alisim_model,
        }

    def load_all_models(self):
        models = {}
        for path in sorted(self.iqtree_dir.glob("*.iqtree")):
            row = self.parse_iqtree(path)
            models[row["pfam_id"]] = row
        return models

    def write_species_model_csv(self, species, pfams, models):
        out_csv = self.species_model_dir / f"{species}_pfam_models.csv"
        with open(out_csv, "w", newline="") as out:
            writer = csv.DictWriter(
                out,
                fieldnames=["pfam_id", "model", "I", "alpha", "rates", "alisim_model"],
            )
            writer.writeheader()
            for pfam_id in sorted(pfams):
                model_row = models.get(pfam_id)
                if model_row is None:
                    writer.writerow({
                        "pfam_id": pfam_id,
                        "model": "",
                        "I": "",
                        "alpha": "",
                        "rates": "",
                        "alisim_model": "",
                    })
                else:
                    writer.writerow({
                        "pfam_id": pfam_id,
                        "model": model_row["model"],
                        "I": model_row["I"],
                        "alpha": model_row["alpha"],
                        "rates": model_row["rates"],
                        "alisim_model": model_row["alisim_model"],
                    })
        return out_csv

    def write_all_species_model_csvs(self, species_pfams, models):
        written = []
        for species in sorted(species_pfams):
            out_csv = self.write_species_model_csv(
                species=species,
                pfams=species_pfams[species],
                models=models,
            )
            written.append(out_csv)
        return written

    def run(self):
        self.make_dirs()
        print(f"PFAM scan folder: {self.pfam_scan_folder}")
        print(f"Results folder: {self.results_dir}")
        print(f"Species CSV folder: {self.species_model_dir}")
        print(f"IQ-TREE threads: {self.threads}")

        species_pfams, all_pfams = self.collect_species_pfams()
        print(f"\nFound {len(species_pfams)} PFAM scan files")
        print(f"Found {len(all_pfams)} unique PFAM domains total")

        for i, pfam_id in enumerate(all_pfams, start=1):
            print(f"\n[{i}/{len(all_pfams)}] Processing {pfam_id}")
            try:
                self.process_one_pfam(pfam_id)
            except subprocess.CalledProcessError as err:
                print(f"WARNING: command failed for {pfam_id}: {err}")
                continue
            except Exception as err:
                print(f"WARNING: unexpected failure for {pfam_id}: {err}")
                continue

        models = self.load_all_models()
        if all_pfams and not any(p in models for p in all_pfams):
            # without any domain model the species cannot be simulated, so stop here
            # rather than write a model table with every model empty
            print("\nERROR: no substitution model could be estimated for any of the "
                  f"{len(all_pfams)} PFAM domains (see the warnings above).")
            sys.exit(1)
        written_csvs = self.write_all_species_model_csvs(species_pfams, models)
        print("\nDone.")
        print(f"Wrote {len(written_csvs)} per-species model CSVs to:")
        print(self.species_model_dir)


def main():
    args = parse_args()
    threads = safe_int_threads(args.threads)
    extractor = PfamModelExtractor(
        pfam_scan_folder=args.pfam_scan_folder,
        results_dir=args.results_dir,
        threads=threads,
    )
    extractor.run()


if __name__ == "__main__":
    main()
