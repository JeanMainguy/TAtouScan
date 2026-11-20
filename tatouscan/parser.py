from pathlib import Path
from typing import DefaultDict, Dict, List
from collections import defaultdict
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
            key, value = attribute.strip().split("=")
            attributes_dict[key.upper()] = value
        except ValueError:
            pass  # we assume that it is a strange, but useless field for our analysis

    return attributes_dict


def get_cdss_from_gff_file(gff_file: Path):
    """Parse a GFF file and return a list of GFF entries."""

    contig_id_to_cds: DefaultDict[str, List[Cds]] = defaultdict(list)
    with read_file(gff_file) as gff_fh:

        for line_number, line in enumerate(gff_fh, start=1):

            if line.startswith("##FASTA"):
                # Sequences appended after the annotations (e.g. Prokka, Bakta)
                break

            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"Line {line_number} of {gff_file} has {len(fields)} tab-separated "
                    f"columns instead of 9: {line.rstrip()!r}"
                )

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
            ) = fields

            if feature in ["CDS", "region"]:

                if feature == "region":
                    # retrieve info on the contig
                    pass

                elif feature == "CDS":

                    attributes = parse_gff_attributes(attributes)

                    gene_id = attributes.get("ID")
                    protein_id = attributes.get("PROTEIN_ID")
                    locus_tag = attributes.get("LOCUS_TAG")
                    name = attributes.get("NAME")
                    product = attributes.get("PRODUCT")

                    if gene_id is None:
                        raise ValueError(
                            f"Missing ID attribute in CDS feature: {line}, of file {gff_file}"
                        )

                    if strand not in ("+", "-"):
                        raise ValueError(
                            f"CDS '{gene_id}' at line {line_number} of {gff_file} has "
                            f"strand '{strand}': expected '+' or '-'."
                        )

                    cds = Cds(
                        id=gene_id,
                        contig_id=contig_id,
                        coordinates=[(int(start), int(stop))],
                        strand=Strand(strand),
                        # Some tools write "." instead of 0/1/2 as CDS phase
                        frame=None if frame == "." else Frame(int(frame)),
                        product=product,
                        locus_tag=locus_tag,
                        name=name,
                        protein_id=protein_id,
                    )
                    contig_id_to_cds[contig_id].append(cds)

    for contig, cds_list in contig_id_to_cds.items():
        yield contig, cds_list
