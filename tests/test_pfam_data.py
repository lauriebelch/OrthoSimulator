'''
Checks that the precomputed PFAM data for the twelve model species is present
and readable.

The data may be present as plain files or only as the compressed archive
PFAM/archives/<species>.tar.xz, which OrthoSimulator unpacks on first use.
The tests below unpack into a temporary copy, so they work either way.

For each species, OrthoSimulator needs three files:
  PFAM/pfam_genomes/<species>.fa
  PFAM/pfamscan_results/<species>_pfam.txt
  PFAM/pfam_results/csv_files/species_models_fixed/<species>_pfam_models.csv
'''

import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import csv
import io
import tarfile
import tempfile

from scripts.PFAM import pfam_archive

PFAM_ROOT = os.path.join(REPO_ROOT, "PFAM")

# the twelve species named in the paper
SPECIES = [
    "Homo_sapiens",
    "Arabidopsis_thaliana",
    "Saccharomyces_cerevisiae",
    "Escherichia_coli",
    "Mnemiopsis_leidyi",
    "Danio_rerio",
    "Gallus_gallus",
    "Mus_musculus",
    "Chlamydomonas_reinhardtii",
    "Drosophila_melanogaster",
    "Triticum_aestivum",
    "Haloferax_volcanii",
]


_READY_ROOT = None


def ready_root():
    """
    A PFAM folder that has every species unpacked. If the plain files are all
    in the repo this is the repo folder itself; otherwise the archives are
    unpacked once into a temporary folder.
    """
    global _READY_ROOT
    if _READY_ROOT is None:
        plain = all(os.path.isfile(os.path.join(PFAM_ROOT, rel))
                    for s in SPECIES for rel in pfam_archive.species_files(s))
        if plain:
            _READY_ROOT = PFAM_ROOT
        else:
            tmp = tempfile.mkdtemp(prefix="pfam_test_")
            os.makedirs(os.path.join(tmp, "archives"))
            for s in SPECIES:
                src = pfam_archive.archive_path(PFAM_ROOT, s)
                assert os.path.isfile(src), f"neither plain files nor an archive for {s}"
                dst = pfam_archive.archive_path(tmp, s)
                os.symlink(src, dst)
                pfam_archive.unpack(tmp, s)
            _READY_ROOT = tmp
    return _READY_ROOT


def genome_path(species):
    return os.path.join(ready_root(), "pfam_genomes", f"{species}.fa")


def scan_path(species):
    return os.path.join(ready_root(), "pfamscan_results", f"{species}_pfam.txt")


def models_path(species):
    return os.path.join(
        ready_root(), "pfam_results", "csv_files", "species_models_fixed",
        f"{species}_pfam_models.csv",
    )


def test_all_twelve_species_have_all_three_files():
    missing = []
    for species in SPECIES:
        for path in (genome_path(species), scan_path(species), models_path(species)):
            if not os.path.isfile(path):
                missing.append(os.path.relpath(path, REPO_ROOT))
    assert missing == [], "missing PFAM files: " + ", ".join(missing)


def test_archives_exist_for_all_twelve_species():
    for species in SPECIES:
        assert os.path.isfile(pfam_archive.archive_path(PFAM_ROOT, species)), f"no archive for {species}"


def test_genome_files_are_protein_fasta():
    for species in SPECIES:
        with open(genome_path(species)) as f:
            first = f.readline()
        assert first.startswith(">"), f"{species}: genome file does not start with a FASTA header"


def test_pfam_scan_files_contain_domain_hits():
    # pfam_scan.pl output: comment lines start with '#', hit rows have the
    # hit type ('Domain', 'Family', ...) in column 8
    for species in SPECIES:
        found_domain = False
        with open(scan_path(species)) as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 8 and parts[7] == "Domain":
                    found_domain = True
                    break
        assert found_domain, f"{species}: no Domain hits in the pfam_scan file"


def test_model_tables_have_alisim_models():
    for species in SPECIES:
        with open(models_path(species), newline="") as f:
            reader = csv.DictReader(f)
            assert reader.fieldnames is not None and "alisim_model" in reader.fieldnames, \
                f"{species}: no 'alisim_model' column"
            assert any(row.get("alisim_model") for row in reader), \
                f"{species}: the 'alisim_model' column is empty"


def test_no_private_paths_in_pfam_files():
    # comment headers in pfam_scan output can contain the path of whoever ran it
    for species in SPECIES:
        with open(scan_path(species)) as f:
            for line in f:
                if not line.startswith("#"):
                    break
                assert "/local/home" not in line and "/Users/" not in line, \
                    f"{species}: private path in pfam_scan header"


def test_unpack_round_trip_and_does_not_overwrite():
    with tempfile.TemporaryDirectory() as root:
        contents = {}
        for rel in pfam_archive.species_files("Toy_species"):
            os.makedirs(os.path.dirname(os.path.join(root, rel)), exist_ok=True)
            contents[rel] = ("data for " + rel + "\n") * 50
            with open(os.path.join(root, rel), "w") as f:
                f.write(contents[rel])
        pfam_archive.pack(root, "Toy_species")
        for rel in contents:
            os.remove(os.path.join(root, rel))
        assert pfam_archive.unpack(root, "Toy_species") == 3
        for rel, text in contents.items():
            with open(os.path.join(root, rel)) as f:
                assert f.read() == text
        # a second call finds everything in place and writes nothing
        assert pfam_archive.unpack(root, "Toy_species") == 0
        # a species with no archive is not an error
        assert pfam_archive.unpack(root, "No_such_species") == 0


def test_unpack_rejects_files_outside_the_expected_names():
    with tempfile.TemporaryDirectory() as root:
        os.makedirs(os.path.join(root, "archives"))
        with tarfile.open(pfam_archive.archive_path(root, "Bad"), "w:xz") as tar:
            data = b"x"
            info = tarfile.TarInfo("../escape.txt")
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
        try:
            pfam_archive.unpack(root, "Bad")
        except ValueError:
            pass
        else:
            raise AssertionError("an unexpected archive entry was accepted")
        assert not os.path.exists(os.path.join(root, "..", "escape.txt"))


def test_unpack_is_safe_when_many_processes_unpack_the_same_species_at_once():
    '''The GA evaluates solutions in parallel, and on a fresh copy each one unpacks the
    starting species the first time it is used. On the server (8 Nov) several processes
    wrote the same temporary file and all but one failed to rename it.'''
    import multiprocessing as mp
    species = "Saccharomyces_cerevisiae"
    with tempfile.TemporaryDirectory() as tmp:
        os.makedirs(os.path.dirname(pfam_archive.archive_path(tmp, species)), exist_ok=True)
        with open(pfam_archive.archive_path(PFAM_ROOT, species), "rb") as src, \
             open(pfam_archive.archive_path(tmp, species), "wb") as dst:
            dst.write(src.read())
        with mp.get_context("fork").Pool(8) as pool:
            try:
                pool.starmap(pfam_archive.unpack, [(tmp, species)] * 16)
            except Exception as error:
                raise AssertionError(f"parallel unpack failed: {error!r}")
        # every file complete and identical to a single unpack, and no temporary files left
        with tempfile.TemporaryDirectory() as ref:
            os.makedirs(os.path.dirname(pfam_archive.archive_path(ref, species)), exist_ok=True)
            with open(pfam_archive.archive_path(tmp, species), "rb") as src, \
                 open(pfam_archive.archive_path(ref, species), "wb") as dst:
                dst.write(src.read())
            pfam_archive.unpack(ref, species)
            for rel in pfam_archive.species_files(species):
                with open(os.path.join(tmp, rel), "rb") as a, open(os.path.join(ref, rel), "rb") as b:
                    assert a.read() == b.read(), f"{rel} differs after a parallel unpack"
        leftovers = [f for d, _, fs in os.walk(tmp) for f in fs if f.endswith(".part")]
        assert leftovers == [], leftovers


def test_pfam_model_extractor_uses_the_bundled_iqtree():
    '''--PFAM estimates a model per domain with IQ-TREE. It called "iqtree3" from the PATH,
    but IQ-TREE ships in bin/ and is not on the PATH, so every model failed silently.'''
    from scripts.PFAM import PfamModelExtractor as pme
    assert pme.iqtree_executable() == os.path.join(REPO_ROOT, "bin", "iqtree3")


def test_pfam_model_extractor_stops_when_no_domain_model_could_be_made():
    '''With no model at all the new species cannot be simulated (the simulation crashed
    later with an IndexError), so --PFAM must stop with an error instead.'''
    from scripts.PFAM import PfamModelExtractor as pme
    with tempfile.TemporaryDirectory() as tmp:
        scans = os.path.join(tmp, "scans")
        os.makedirs(scans)
        with open(os.path.join(scans, "Toy_pfam.txt"), "w") as f:
            f.write("# header\n\nseq1  1 50 1 50 PF00069.30 Pkinase Domain 1 50 50 40.0 1e-9 1 No_clan\n")
        ex = pme.PfamModelExtractor(scans, os.path.join(tmp, "results"), 1)
        ex.process_one_pfam = lambda pfam_id: None   # every model fails, as when IQ-TREE is missing
        try:
            ex.run()
        except SystemExit as stop:
            assert stop.code == 1
        else:
            raise AssertionError("run() should stop when no domain model could be made")


def test_simulation_says_clearly_when_a_species_has_no_domain_models():
    from types import SimpleNamespace
    import scripts.OrthogroupSimulation.orthogroup_simulation_utils as og
    with tempfile.TemporaryDirectory() as tmp:
        table = os.path.join(tmp, "Toy_pfam_models.csv")
        with open(table, "w") as f:
            f.write("file,model,I,alpha,rates,alisim_model\ngene1,,,,,\ngene2,,,,,\n")
        try:
            og.SelectRandomDomainModel(SimpleNamespace(d=table))
        except ValueError as error:
            assert "Run --PFAM again" in str(error)
        else:
            raise AssertionError("an empty model table should give a clear error")
        with open(table, "a") as f:
            f.write("gene3,JTT+G4,,1.0,,JTT+G4{1.0}\n")
        assert og.SelectRandomDomainModel(SimpleNamespace(d=table)) == "JTT+G4{1.0}"
