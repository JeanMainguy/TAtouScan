from tatouscan.models import Gene, Contig, Cds


def test_gene_initialization():
    gene = Gene(name="gene1", coordinates=[(1, 4)])

    assert gene.name == "gene1"
    assert gene.coordinates == [(1, 4)]
    assert gene.contig is None
    assert gene.start == 1
    assert gene.stop == 4


def test_gene_initialization_with_join_coordinates():
    gene = Gene(name="gene1", coordinates=[(1, 4), (5, 8)])

    assert gene.name == "gene1"
    assert gene.coordinates == [(1, 4), (5, 8)]
    assert gene.contig is None
    assert gene.start == 1
    assert gene.stop == 8


def test_cds_initialization():
    gene = Cds(name="cds1", protein_sequence="MTDLA", coordinates=[(1, 25)])

    assert gene.name == "cds1"
    assert gene.protein_sequence == "MTDLA"
    assert gene.coordinates == [(1, 25)]
    assert gene.contig is None
    assert gene.start == 1
    assert gene.stop == 25


def test_contig_initialization():
    gene1 = Gene(name="gene1", coordinates=[(1, 4), (5, 8)])
    gene2 = Gene(name="gene2", coordinates=[(8, 15)])
    contig = Contig(name="contig1", genes=[gene1, gene2])

    assert contig.name == "contig1"
    assert len(contig.genes) == 2
    assert contig.genes[0].name == "gene1"
    assert contig.genes[1].name == "gene2"
    assert gene1.contig == contig
    assert gene2.contig == contig
