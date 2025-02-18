from enum import Enum
from typing import List, Optional, Tuple

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

    @property
    def protein_sequence(self) -> str:
        """
        stop coordinate of the feature.
        """
        if self.digit_sequence is None:
            raise ValueError(f"CDS {self.id} has no protein sequence.")
        return self.digit_sequence.textize().sequence

    def __repr__(self):
        return f"CDS(name={self.id}, coordinates={self.coordinates}, strand={self.strand}, protein_id={self.protein_id})"

class TaHit:

    def __init__(self, protein_id: str, ta_name: str, score: float, evalue: float):
        self.protein_id = protein_id
        self.ta_name = ta_name
        self.score = score
        self.evalue = evalue
