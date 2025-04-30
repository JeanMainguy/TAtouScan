import gzip
from pathlib import Path


def read_file(gff_file: Path):

    if gff_file.suffix == ".gz":
        return gzip.open(gff_file, "rt")
    else:
        return open(gff_file, "r")
