from pathlib import Path
import pytest
from tatouscan.parser import get_cdss_from_gff_file, parse_gff_attributes


@pytest.fixture
def gff_file(tmp_path: Path):

    gff_file = tmp_path / "genomes.gff"

    with open(gff_file, "w") as fh:
        fh.write(
            "##gff-version 3\n"
            "##sequence-region contig1 1 300\n"
            "##sequence-region contig2 1 500\n"
            "contig1	Prodigal:2.6	CDS	1	25	.	+	0	ID=cds1\n"
            "contig1	Prodigal:2.6	CDS	50	100	.	+	0	ID=cds2\n"
            "contig2	Prodigal:2.6	CDS	1	25	.	+	0	ID=cds3\n"
            "contig2	Prodigal:2.6	CDS	50	100	.	+	0	ID=cds4\n"
        )
    return gff_file


@pytest.fixture
def faa_file(tmp_path: Path):

    gff_file = tmp_path / "genomes.faa"

    with open(gff_file, "w") as fh:
        fh.write(
            ">cds1\n"
            "TGPYMMNA\n"
            ">cds2\n"
            "TGPYMMNA\n"
            ">cds3\n"
            "TGPYMMNA\n"
            ">cds4\n"
            "TGPYMMNA\n"
        )
    return gff_file


def test_parse_gff_file(gff_file: Path):

    contig_name_and_cdss = get_cdss_from_gff_file(gff_file=gff_file)

    contig_name_and_cdss = list(contig_name_and_cdss)

    assert len(contig_name_and_cdss) == 2
    assert contig_name_and_cdss[0][0] == "contig1"
    assert len(contig_name_and_cdss[0][1]) == 2
    assert contig_name_and_cdss[1][0] == "contig2"
    assert len(contig_name_and_cdss[1][1]) == 2


def test_parse_gff_file_with_missing_id(gff_file: Path):

    with open(gff_file, "a") as fh:
        fh.write("contig2	Prodigal:2.6	CDS	50	100	.	+	0	Name=cds5\n")

    with pytest.raises(ValueError):
        list(get_cdss_from_gff_file(gff_file=gff_file))


def test_parse_gff_file_stops_at_fasta_section(gff_file: Path):
    # Prokka and Bakta append the genome sequence after a ##FASTA line
    with open(gff_file, "a") as fh:
        fh.write("##FASTA\n>contig1\nACGTACGT\n>contig2\nACGT\n")

    contig_name_and_cdss = list(get_cdss_from_gff_file(gff_file=gff_file))

    assert [contig for contig, _ in contig_name_and_cdss] == ["contig1", "contig2"]
    assert sum(len(cdss) for _, cdss in contig_name_and_cdss) == 4


def test_parse_gff_file_skips_blank_lines(gff_file: Path):
    with open(gff_file, "a") as fh:
        fh.write("\n   \ncontig2	Prodigal:2.6	CDS	200	300	.	+	0	ID=cds5\n\n")

    contig_name_and_cdss = dict(get_cdss_from_gff_file(gff_file=gff_file))

    assert [cds.id for cds in contig_name_and_cdss["contig2"]] == [
        "cds3",
        "cds4",
        "cds5",
    ]


def test_parse_gff_file_with_undefined_phase(gff_file: Path):
    with open(gff_file, "a") as fh:
        fh.write("contig2	Prodigal:2.6	CDS	200	300	.	+	.	ID=cds5\n")

    contig_name_and_cdss = dict(get_cdss_from_gff_file(gff_file=gff_file))

    cds5 = contig_name_and_cdss["contig2"][-1]
    assert cds5.id == "cds5"
    assert cds5.frame is None


def test_parse_gff_file_with_wrong_column_count(gff_file: Path):
    with open(gff_file, "a") as fh:
        fh.write("contig2	Prodigal:2.6	CDS	200	300	.	+	0\n")

    with pytest.raises(ValueError, match="instead of 9"):
        list(get_cdss_from_gff_file(gff_file=gff_file))


def test_parse_gff_file_with_undefined_strand(gff_file: Path):
    with open(gff_file, "a") as fh:
        fh.write("contig2	Prodigal:2.6	CDS	200	300	.	.	0	ID=cds5\n")

    with pytest.raises(ValueError, match="expected '\\+' or '-'"):
        list(get_cdss_from_gff_file(gff_file=gff_file))


def test_parse_gff_attributes_decodes_percent_encoding():
    attributes = parse_gff_attributes(
        "ID=cds1;product=toxin%2C RelE%2FParE family%3B type II;Note=50%25 identity"
    )

    assert attributes["ID"] == "cds1"
    assert attributes["PRODUCT"] == "toxin, RelE/ParE family; type II"
    assert attributes["NOTE"] == "50% identity"


def test_parse_gff_attributes_keeps_unescaped_equal_sign():
    attributes = parse_gff_attributes("ID=cds1;product=toxin a=b")

    assert attributes["PRODUCT"] == "toxin a=b"


def test_parse_gff_file_decodes_product(gff_file: Path):
    with open(gff_file, "a") as fh:
        fh.write(
            "contig2	Prodigal:2.6	CDS	200	300	.	+	0	ID=cds5;product=toxin%2C RelE family\n"
        )

    contig_name_and_cdss = dict(get_cdss_from_gff_file(gff_file=gff_file))

    assert contig_name_and_cdss["contig2"][-1].product == "toxin, RelE family"
