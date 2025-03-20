from pathlib import Path
from typing import DefaultDict, Dict, List, Any
from collections import defaultdict
from pyfastx import Fasta  # type: ignore
from pyhmmer.easel import Alphabet, TextSequence
from tatouscan.utils import read_file
from tatouscan.models import Cds, Strand, Frame


def parse_gff_attributes(attribute_str: str) -> Dict[str, str]:
    """Parses the gff attribute's line and outputs the attributes_get in a dict structure.

    :param attribute_str: The attribute line from the GFF file.
    :return: A dictionary of the attributes.
    """
    attributes = [field for field in attribute_str.strip().split(";") if len(field) > 0]

    attributes_dict: Dict[str, str] = {}
    for attribute in attributes:
        try:
            (key, value) = attribute.strip().split("=")
            attributes_dict[key.upper()] = value
        except ValueError:
            pass  # we assume that it is a strange, but useless field for our analysis

    return attributes_dict


def parse_gff_file(gff_file: Path):
    """Parse a GFF file and return a list of GFF entries."""

    contig_id_to_cds: DefaultDict[str, List[Cds]] = defaultdict(list)
    with read_file(gff_file) as gff_fh:

        for line in gff_fh:

            if line.startswith("#"):  # Skip comments
                continue

            (
                contig_id,
                _source,
                feature,
                start,
                stop,
                _score,
                strand,
                frame,
                attributes,
            ) = line.split("\t")

            if feature in ["CDS", "region"]:

                if feature == "region":
                    # retrieve info on the contig
                    pass

                elif feature == "CDS":

                    attributes = parse_gff_attributes(attributes)

                    gene_id = attributes.get("ID")

                    if gene_id is None:
                        raise ValueError(
                            f"Missing ID attribute in CDS feature: {line}, of file {gff_file}"
                        )

                    cds = Cds(
                        id=gene_id,
                        contig_id=contig_id,
                        coordinates=[(int(start), int(stop))],
                        strand=Strand(strand),
                        frame=Frame(int(frame)),
                    )

                    contig_id_to_cds[contig_id].append(cds)

    for contig, cds_list in contig_id_to_cds.items():
        yield contig, cds_list


def get_cds_from_gff_and_faa_files(gff_file: Path, faa_file: Path):

    fa: Any = Fasta(faa_file.as_posix(), build_index=True)

    for contig_name, cds_list in parse_gff_file(gff_file):
        for cds in cds_list:
            seq = fa[cds.id]
            sequence = seq.seq.replace("*", "").encode("UTF-8")
            digit_seq = TextSequence(sequence).digitize(Alphabet.amino())
            cds.digit_sequence = digit_seq

        yield contig_name, cds_list
