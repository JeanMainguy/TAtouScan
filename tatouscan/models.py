from enum import Enum
from typing import List, Optional, Tuple, Set

from pyhmmer.easel import DigitalSequence


class Strand(Enum):
    POSITIVE = "+"
    NEGATIVE = "-"


class Frame(Enum):
    ZERO = 0
    ONE = 1
    TWO = 2


class Contig:

    def __init__(
        self,
        id: str,
        genes: List["Gene"] = [],
        length: Optional[int] = None,
    ):
        self.id = id
        self.genes = genes or []
        self.length = length

    def __repr__(self):
        return f"Contig(id={self.id}, num_genes={len(self.genes)})"


class Gene:

    def __init__(
        self,
        id: str,
        contig_id: str,
        coordinates: list[Tuple[int, int]],
        strand: Strand,
        product: Optional[str] = None,
        locus_tag: Optional[str] = None,
        name: Optional[str] = None,
    ):
        self.id = id
        self.coordinates = coordinates
        self.strand = strand
        self.contig_id = contig_id
        self.name = name
        self.product = product
        self.locus_tag = locus_tag

    def __repr__(self):
        return f"Gene(id={self.id}, coordinates={self.coordinates})"

    @property
    def start(self) -> int:
        """
        start coordinate of the feature.
        """
        return self.coordinates[0][0]

    @property
    def stop(self) -> int:
        """
        stop coordinate of the feature.
        """
        return self.coordinates[-1][-1]

    def distance_from(self, other_gene: "Gene"):
        # TODO take into account circularity
        assert self.contig_id == other_gene.contig_id

        if self.start <= other_gene.start:
            return other_gene.start - self.stop

        elif self.start > other_gene.start:
            return self.start - other_gene.stop
        else:
            # self.start == other_gene.start
            raise ValueError()


class Cds(Gene):

    def __init__(
        self,
        id: str,
        contig_id: str,
        coordinates: list[Tuple[int, int]],
        strand: Strand,
        frame: Frame,
        product: Optional[str] = None,
        locus_tag: Optional[str] = None,
        name: Optional[str] = None,
        protein_id: Optional[str] = None,
        digit_sequence: Optional[DigitalSequence] = None,
    ):
        super().__init__(
            id=id,
            contig_id=contig_id,
            coordinates=coordinates,
            strand=strand,
            product=product,
            locus_tag=locus_tag,
            name=name,
        )

        self.frame = frame
        self.protein_id = protein_id
        self.digit_sequence = digit_sequence
        self.ta_hits: List[TaHit] = []
        self.neighbor_genes: List["Cds"] = []
        self.ta_cluster: GeneCluster | None = None

    def __repr__(self):
        return f"CDS(name={self.id}, coordinates={self.coordinates}, strand={self.strand}, protein_id={self.protein_id})"

    def add_neigbor_gene(self, gene: "Cds"):
        self.neighbor_genes.append(gene)
        gene.neighbor_genes.append(self)

    def update_ta_cluster(self, gene: "Cds"):

        if self.ta_cluster is None and gene.ta_cluster is None:

            ta_cluster: GeneCluster = GeneCluster({self, gene})

            self.ta_cluster = ta_cluster
            gene.ta_cluster = ta_cluster

        elif self.ta_cluster and gene.ta_cluster is None:

            self.ta_cluster.add(gene)
            gene.ta_cluster = self.ta_cluster

        elif gene.ta_cluster and self.ta_cluster is None:
            gene.ta_cluster.add(self)
            self.ta_cluster = gene.ta_cluster

        elif gene.ta_cluster and self.ta_cluster:

            self.ta_cluster.combine(gene.ta_cluster)


class GeneCluster:
    counter: int = 1

    def __init__(self, genes: Set[Cds]):
        self.genes = genes

        self.id = GeneCluster.counter
        GeneCluster.counter += 1

    def add(self, gene: Cds):
        self.genes.add(gene)

    def combine(self, other_cluster: "GeneCluster"):

        self.genes |= other_cluster.genes

        for gene in self.genes:
            gene.ta_cluster = self

    def __repr__(self):
        return f"GeneCluster(gene_count={len(self.genes)}: {[gene for gene in self.genes]})"


class TaHit:
    def __init__(self, protein_id: str, ta_name: str, score: float, evalue: float):
        self.protein_id = protein_id
        self.ta_name = ta_name
        self.score = score
        self.evalue = evalue

    def __repr__(self):
        return f"TaHit(protein_id={self.protein_id}, ta_name={self.ta_name}, score={self.score}, evalue={self.evalue})"
