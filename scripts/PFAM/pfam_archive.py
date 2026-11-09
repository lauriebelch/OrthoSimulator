#!/usr/bin/env python3
'''
Pack and unpack the precomputed PFAM data, one compressed archive per species.

Each species needs three files (paths are relative to the PFAM folder):
  pfam_genomes/<species>.fa
  pfamscan_results/<species>_pfam.txt
  pfam_results/csv_files/species_models_fixed/<species>_pfam_models.csv

They are stored together in PFAM/archives/<species>.tar.xz. OrthoSimulator
unpacks a species the first time it is used, so users do not need to do
anything. To unpack everything in advance:

  python3 scripts/PFAM/pfam_archive.py unpack            # all species
  python3 scripts/PFAM/pfam_archive.py unpack Homo_sapiens

To rebuild the archives from the plain files (maintainers only):

  python3 scripts/PFAM/pfam_archive.py pack
'''

import argparse
import glob
import os
import sys
import tarfile
import uuid

ARCHIVE_DIR = "archives"


def species_files(species):
    return [
        os.path.join("pfam_genomes", f"{species}.fa"),
        os.path.join("pfamscan_results", f"{species}_pfam.txt"),
        os.path.join("pfam_results", "csv_files", "species_models_fixed", f"{species}_pfam_models.csv"),
    ]


def archive_path(pfam_root, species):
    return os.path.join(pfam_root, ARCHIVE_DIR, f"{species}.tar.xz")


def archived_species(pfam_root):
    paths = glob.glob(os.path.join(pfam_root, ARCHIVE_DIR, "*.tar.xz"))
    return sorted(os.path.basename(p)[:-len(".tar.xz")] for p in paths)


def pack(pfam_root, species):
    os.makedirs(os.path.join(pfam_root, ARCHIVE_DIR), exist_ok=True)

    def clean(info):
        # same archive bytes every time the same files are packed
        info.mtime = 0
        info.uid = info.gid = 0
        info.uname = info.gname = ""
        return info

    with tarfile.open(archive_path(pfam_root, species), "w:xz", preset=6) as tar:
        for rel in species_files(species):
            tar.add(os.path.join(pfam_root, rel), arcname=rel, filter=clean)


def unpack(pfam_root, species):
    '''Unpack the files for one species that are not already there.
    Returns the number of files written (0 if there is no archive or nothing was missing).'''
    path = archive_path(pfam_root, species)
    written = 0
    if not os.path.isfile(path):
        return written
    wanted = species_files(species)
    with tarfile.open(path, "r:xz") as tar:
        for member in tar.getmembers():
            # only plain files with the expected names, nothing outside pfam_root
            if member.name not in wanted or not member.isfile():
                raise ValueError(f"unexpected entry in {path}: {member.name}")
            target = os.path.join(pfam_root, member.name)
            if os.path.exists(target):
                continue
            os.makedirs(os.path.dirname(target), exist_ok=True)
            # a temporary name of its own: several processes (for example the GA's
            # parallel fitness evaluations) can unpack the same species at once. Each
            # writes a complete copy and renames it into place; the rename is atomic,
            # so the target is never partly written and the last identical copy wins.
            tmp = f"{target}.{os.getpid()}.{uuid.uuid4().hex}.part"
            try:
                with tar.extractfile(member) as src, open(tmp, "wb") as dst:
                    while True:
                        block = src.read(1 << 20)
                        if not block:
                            break
                        dst.write(block)
                os.replace(tmp, target)
            finally:
                if os.path.exists(tmp):
                    os.remove(tmp)
            written += 1
    return written


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("command", choices=["pack", "unpack"])
    p.add_argument("species", nargs="*", help="default: all species")
    p.add_argument("--pfam-root", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "PFAM"))
    args = p.parse_args()
    root = os.path.abspath(args.pfam_root)

    if args.command == "pack":
        species = args.species or sorted(
            os.path.basename(f)[:-3] for f in glob.glob(os.path.join(root, "pfam_genomes", "*.fa")))
        for s in species:
            pack(root, s)
            print(f"packed {s}: {os.path.getsize(archive_path(root, s)) / 1e6:.1f} MB")
    else:
        species = args.species or archived_species(root)
        for s in species:
            if not os.path.isfile(archive_path(root, s)):
                sys.exit(f"no archive for {s} in {os.path.join(root, ARCHIVE_DIR)}")
            print(f"{s}: {unpack(root, s)} file(s) unpacked")


if __name__ == "__main__":
    main()
