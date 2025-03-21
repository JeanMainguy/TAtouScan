from pathlib import Path
import pytest
from tatouscan.parser import parse_gff_file, get_cds_from_gff_and_faa_files


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

    contig_name_and_cdss = parse_gff_file(gff_file=gff_file)

    contig_name_and_cdss = list(contig_name_and_cdss)

    assert len(contig_name_and_cdss) == 2
    assert contig_name_and_cdss[0][0] == "contig1"
    assert len(contig_name_and_cdss[0][1]) == 2
    assert contig_name_and_cdss[1][0] == "contig2"
    assert len(contig_name_and_cdss[1][1]) == 2


def test_parse_gff_file_with_missing_id(gff_file: Path):

    with open(gff_file, "a") as fh:
        fh.write("contig2	Prodigal:2.6	CDS	50	100	.	+	0\n")

    with pytest.raises(ValueError):
        list(parse_gff_file(gff_file=gff_file))


def test_get_cds_from_gff_and_faa_files(gff_file:Path, faa_file:Path):
    
    contig_name_and_cdss = get_cds_from_gff_and_faa_files(gff_file, faa_file)

    contig_name_and_cdss = list(contig_name_and_cdss)

    assert len(contig_name_and_cdss) == 2

    contig_name, cdss = contig_name_and_cdss[0]
    cds = cdss[0]
    assert contig_name == "contig1"
    assert cds.id == "cds1"
    assert cds.start == 1
    assert cds.stop == 25
    assert cds.protein_sequence == "TGPYMMNA"