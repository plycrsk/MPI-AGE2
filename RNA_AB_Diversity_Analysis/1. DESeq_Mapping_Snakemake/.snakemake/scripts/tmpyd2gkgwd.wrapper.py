
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/home/plycrsk/anaconda3/envs/snakemake/lib/python3.10/site-packages']); import pickle; snakemake = pickle.loads(b'\x80\x04\x95\xc9\x05\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94}\x94(\x8c\x06_names\x94}\x94\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x0f\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x15)}\x94\x8c\x05_name\x94h\x0fsNt\x94bh\x10h\x13h\x15\x85\x94R\x94(h\x15)}\x94h\x19h\x10sNt\x94bub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94\x8c\x16resources/genome.fasta\x94a}\x94(h\x0b}\x94h\r]\x94(h\x0fh\x10eh\x0fh\x13h\x15\x85\x94R\x94(h\x15)}\x94h\x19h\x0fsNt\x94bh\x10h\x13h\x15\x85\x94R\x94(h\x15)}\x94h\x19h\x10sNt\x94bub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94(\x8c\x16nothobranchius_furzeri\x94\x8c\x03dna\x94\x8c\x0cNfu_20140520\x94Khe}\x94(h\x0b}\x94(\x8c\x07species\x94K\x00N\x86\x94\x8c\x08datatype\x94K\x01N\x86\x94\x8c\x05build\x94K\x02N\x86\x94\x8c\x07release\x94K\x03N\x86\x94uh\r]\x94(h\x0fh\x10eh\x0fh\x13h\x15\x85\x94R\x94(h\x15)}\x94h\x19h\x0fsNt\x94bh\x10h\x13h\x15\x85\x94R\x94(h\x15)}\x94h\x19h\x10sNt\x94bh8h3h:h4h<h5h>Khub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94}\x94(h\x0b}\x94h\r]\x94(h\x0fh\x10eh\x0fh\x13h\x15\x85\x94R\x94(h\x15)}\x94h\x19h\x0fsNt\x94bh\x10h\x13h\x15\x85\x94R\x94(h\x15)}\x94h\x19h\x10sNt\x94bub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01\x8c\x04/tmp\x94e}\x94(h\x0b}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06tmpdir\x94K\x02N\x86\x94uh\r]\x94(h\x0fh\x10eh\x0fh\x13h\x15\x85\x94R\x94(h\x15)}\x94h\x19h\x0fsNt\x94bh\x10h\x13h\x15\x85\x94R\x94(h\x15)}\x94h\x19h\x10sNt\x94bh`K\x01hbK\x01hdh]ub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94\x8c\x13logs/get-genome.log\x94a}\x94(h\x0b}\x94h\r]\x94(h\x0fh\x10eh\x0fh\x13h\x15\x85\x94R\x94(h\x15)}\x94h\x19h\x0fsNt\x94bh\x10h\x13h\x15\x85\x94R\x94(h\x15)}\x94h\x19h\x10sNt\x94bub\x8c\x06config\x94}\x94(\x8c\x07samples\x94\x8c\x12config/samples.tsv\x94\x8c\x05units\x94\x8c\x10config/units.tsv\x94\x8c\x03ref\x94}\x94(\x8c\x07species\x94h3\x8c\x07release\x94Kh\x8c\x05build\x94h5u\x8c\x08trimming\x94}\x94\x8c\x08activate\x94\x89s\x8c\nmergeReads\x94}\x94\x8c\x08activate\x94\x89s\x8c\x03pca\x94}\x94(\x8c\x08activate\x94\x88\x8c\x06labels\x94]\x94\x8c\tcondition\x94au\x8c\x07diffexp\x94}\x94(\x8c\tcontrasts\x94}\x94\x8c\x0cold-vs-young\x94]\x94(\x8c\x03old\x94\x8c\x05young\x94es\x8c\x05model\x94\x8c\n~condition\x94u\x8c\x06params\x94}\x94(\x8c\x0bcutadapt-pe\x94\x8c\x00\x94\x8c\x0bcutadapt-se\x94h\xa3\x8c\x04star\x94h\xa3uu\x8c\x04rule\x94\x8c\nget_genome\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8cYhttps://github.com/snakemake/snakemake-wrappers/raw/0.77.0/bio/reference/ensembl-sequence\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = 'https://github.com/snakemake/snakemake-wrappers/raw/0.77.0/bio/reference/ensembl-sequence/wrapper.py';
######## snakemake preamble end #########
__author__ = "Johannes Köster"
__copyright__ = "Copyright 2019, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import subprocess as sp
import sys
from itertools import product
from snakemake.shell import shell

species = snakemake.params.species.lower()
release = int(snakemake.params.release)
build = snakemake.params.build

branch = ""
if release >= 81 and build == "GRCh37":
    # use the special grch37 branch for new releases
    branch = "grch37/"

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

spec = ("{build}" if int(release) > 75 else "{build}.{release}").format(
    build=build, release=release
)

suffixes = ""
datatype = snakemake.params.get("datatype", "")
chromosome = snakemake.params.get("chromosome", "")
if datatype == "dna":
    if chromosome:
        suffixes = ["dna.chromosome.{}.fa.gz".format(chromosome)]
    else:
        suffixes = ["dna.primary_assembly.fa.gz", "dna.toplevel.fa.gz"]
elif datatype == "cdna":
    suffixes = ["cdna.all.fa.gz"]
elif datatype == "cds":
    suffixes = ["cds.all.fa.gz"]
elif datatype == "ncrna":
    suffixes = ["ncrna.fa.gz"]
elif datatype == "pep":
    suffixes = ["pep.all.fa.gz"]
else:
    raise ValueError("invalid datatype, must be one of dna, cdna, cds, ncrna, pep")

if chromosome:
    if not datatype == "dna":
        raise ValueError(
            "invalid datatype, to select a single chromosome the datatype must be dna"
        )

success = False
for suffix in suffixes:
    url = "ftp://ftp.ensembl.org/pub/{branch}release-{release}/fasta/{species}/{datatype}/{species_cap}.{spec}.{suffix}".format(
        release=release,
        species=species,
        datatype=datatype,
        spec=spec.format(build=build, release=release),
        suffix=suffix,
        species_cap=species.capitalize(),
        branch=branch,
    )

    try:
        shell("curl -sSf {url} > /dev/null 2> /dev/null")
    except sp.CalledProcessError:
        continue

    shell("(curl -L {url} | gzip -d > {snakemake.output[0]}) {log}")
    success = True
    break

if not success:
    print(
        "Unable to download requested sequence data from Ensembl. "
        "Did you check that this combination of species, build, and release is actually provided?",
        file=sys.stderr,
    )
    exit(1)
