from tatouscan.models import Cds, Frame, GeneCluster, Strand, TaHit
from tatouscan.scoring import GLOBAL_FAMILY, score_cluster

HMM_DB_INFO = {
    "T1": {"type": "Toxin", "source": "TADB3"},
    "T2": {"type": "Toxin", "source": "TADB3"},
    "AT1": {"type": "Antitoxin", "source": "TADB3"},
}

# T1 and T2 share identical statistics, so both toxins score the same.
_STATS = {
    "n_pairs": 10,
    "toxin_size": (100.0, 10.0),
    "at_size": (80.0, 10.0),
    "intergenic_distance": (0.0, 10.0),
}
REF_STATS = {"T1": dict(_STATS), "T2": dict(_STATS), GLOBAL_FAMILY: dict(_STATS)}


def _cds(id, start, stop, hit_name):
    cds = Cds(
        id=id,
        contig_id="contig1",
        coordinates=[(start, stop)],
        strand=Strand("+"),
        frame=Frame(0),
    )
    cds.ta_hits.append(TaHit(protein_id=id, ta_name=hit_name, score=50, evalue=1e-10))
    return cds


def test_score_cluster_tie_break_is_independent_of_gene_order():
    # Antitoxin sits between two equally sized toxins at equal distance,
    # so both (toxin, antitoxin) pairs get exactly the same score.
    toxin_a = _cds("toxA", 1, 300, "T1")
    antitoxin = _cds("at", 301, 540, "AT1")
    toxin_b = _cds("toxB", 541, 840, "T2")

    results = []
    for order in ([toxin_a, antitoxin, toxin_b], [toxin_b, antitoxin, toxin_a]):
        cluster = GeneCluster(set(order))
        cluster.genes = order  # force a specific iteration order
        results.append(
            score_cluster(cluster, HMM_DB_INFO, REF_STATS, known_pairs=frozenset())
        )

    assert results[0] == results[1]
    # The upstream-most toxin wins the tie.
    assert results[0]["matched_family"] == "T1"
