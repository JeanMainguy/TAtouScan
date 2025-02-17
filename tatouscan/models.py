from typing import List, Optional, Tuple


class Contig:
    def __init__(
        self,
        name: str,
        genes: List["Gene"] = [],
        length: Optional[int] = None,
    ):
        self.name = name
        self.genes = genes or []
        self.length = length

    def __repr__(self):
        return f"Contig(id={self.name}, num_genes={len(self.genes)})"


class Gene:
    def __init__(
        self,
        name: str,
        coordinates: list[Tuple[int, int]],
        contig: Optional[Contig] = None,
    ):
        self.name = name
        self.coordinates = coordinates

    def __repr__(self):
        return (
            f"Gene(name={self.name},"
            f"contig={self.contig.name if self.contig else None})"
        )

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
        name: str,
        coordinates: list[Tuple[int, int]],
        protein_sequence: Optional[str] = None,
        contig: Optional[Contig] = None,
    ):
        super().__init__(name, coordinates, contig=contig)
        self.protein_sequence = protein_sequence

    def __repr__(self):
        return (
            f"Gene(name={self.name}, sequence={self.protein_sequence}, "
            f"contig={self.contig.name if self.contig else None})"
        )
