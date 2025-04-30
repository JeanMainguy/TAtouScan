from tatouscan.models import Gene, Cds, Strand, Frame


def test_gene_initialization():
    gene = Gene(
        id="gene1", coordinates=[(1, 4)], contig_id="contig1", strand=Strand("+")
    )

    assert gene.id == "gene1"
    assert gene.coordinates == [(1, 4)]
    assert gene.start == 1
    assert gene.stop == 4


def test_gene_initialization_with_join_coordinates():
    gene = Gene(
        id="gene1",
        coordinates=[(1, 4), (5, 8)],
        contig_id="contig1",
        strand=Strand("+"),
    )

    assert gene.id == "gene1"
    assert gene.coordinates == [(1, 4), (5, 8)]
    assert gene.start == 1
    assert gene.stop == 8
    assert gene.strand == Strand("+")


def test_cds_initialization():
    gene = Cds(
        id="cds1",
        coordinates=[(1, 25)],
        contig_id="contig1",
        strand=Strand("+"),
        frame=Frame(0),
    )

    assert gene.id == "cds1"
    assert gene.coordinates == [(1, 25)]
    assert gene.start == 1
    assert gene.stop == 25
