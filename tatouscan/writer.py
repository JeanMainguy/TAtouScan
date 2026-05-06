from typing import Any, Generator, Dict, List, Optional

from tatouscan.models import Cds, GeneCluster, TaHit
from pathlib import Path
import csv
import itertools
import logging

from tatouscan.scoring import (
    get_gene_type,
    score_cluster,
    score_pair,
    RefStats,
    KnownPairs,
)

logger = logging.getLogger(__name__)

# Score keys that are only written in --detailed mode
_SCORE_DETAIL_KEYS: frozenset[str] = frozenset(
    {
        "matched_family",
        "n_reference_pairs",
        "toxin_size_z",
        "at_size_z",
        "intergenic_distance_z",
    }
)

# Pairs-file columns omitted in default (non-detailed) output
_PAIRS_DETAIL_KEYS: frozenset[str] = frozenset(
    {
        "toxin_start",
        "toxin_end",
        "antitoxin_start",
        "antitoxin_end",
        "toxin_gene_type",
        "antitoxin_gene_type",
    }
)

_SCORE_NULL_FULL: Dict[str, Any] = {
    "matched_family": None,
    "n_reference_pairs": None,
    "toxin_size_z": None,
    "at_size_z": None,
    "intergenic_distance_z": None,
    "pair_is_known": None,
    "score": None,
}


def get_best_hit(ta_hits: List[TaHit], hmm_db_info: Dict[str, Dict[str, str]]):
    """
    Get the best hit from a list of TaHit objects.
    The best hit is defined as the one with the lowest e-value.
    """
    if not ta_hits:
        return {
            "hmm_name": None,
            "hmm_score": None,
            "hmm_evalue": None,
            "hmm_description": None,
        }

    # Sort hits by e-value and return the first one
    best_hit = min(ta_hits, key=lambda hit: hit.evalue)

    hit_info: Dict[str, Any] = {
        "hmm_name": best_hit.ta_name,
        "hmm_score": best_hit.score,
        "hmm_evalue": best_hit.evalue,
        "hmm_description": None,
    }

    # Merge in additional metadata from hmm_db_info using ta_name as key
    if best_hit.ta_name in hmm_db_info:
        hit_info["hmm_description"] = hmm_db_info[best_hit.ta_name][
            "supplementary_info"
        ]

    return hit_info


def _summarise_ta_hits_simple(
    ta_hits: List[TaHit], hmm_db_info: Dict[str, Dict[str, str]]
) -> Dict[str, Any]:
    """Return gene type and the single best HMM hit across all database sources."""
    best_hit = min(ta_hits, key=lambda hit: hit.evalue)
    info = hmm_db_info.get(best_hit.ta_name, {})
    return {
        "gene_type": info.get("type"),
        "hmm_name": best_hit.ta_name,
        "hmm_score": best_hit.score,
        "hmm_evalue": best_hit.evalue,
        "hmm_source": info.get("source"),
        "hmm_description": info.get("supplementary_info"),
    }


def summarise_ta_hits(ta_hits: List[TaHit], hmm_db_info: Dict[str, Dict[str, str]]):

    best_hit = min(ta_hits, key=lambda hit: hit.evalue)
    gene_type = hmm_db_info[best_hit.ta_name]["type"]

    tadb_hits = [
        ta_hit
        for ta_hit in ta_hits
        if ta_hit.ta_name in hmm_db_info
        and hmm_db_info[ta_hit.ta_name]["source"] == "TADB3"
    ]

    best_tadb3_hit = get_best_hit(tadb_hits, hmm_db_info)

    tasmania_hits = [
        ta_hit
        for ta_hit in ta_hits
        if ta_hit.ta_name in hmm_db_info
        and hmm_db_info[ta_hit.ta_name]["source"] == "TASmania"
    ]

    best_tasmania_hit = get_best_hit(tasmania_hits, hmm_db_info)

    other_hits = [
        ta_hit
        for ta_hit in ta_hits
        if ta_hit.ta_name in hmm_db_info
        and hmm_db_info[ta_hit.ta_name]["source"] not in ["TADB3", "TASmania"]
    ]

    best_other_hit = get_best_hit(other_hits, hmm_db_info)

    hit_info: Dict[str, Any] = {
        "gene_type": gene_type,
    }

    hit_info.update({f"TASmania_{k}": v for k, v in best_tasmania_hit.items()})
    hit_info.update({f"Other_{k}": v for k, v in best_other_hit.items()})
    hit_info.update({f"TADB3_{k}": v for k, v in best_tadb3_hit.items()})

    return hit_info


def write_gene_with_ta_annotation(
    contig_name_and_cdss_with_ta_hit: Generator[tuple[str, list[Cds]], Any, None],
    hmm_db_info: Dict[str, Dict[str, str]],
    output_file: Path,
    ref_stats: Optional[RefStats] = None,
    known_pairs: Optional[KnownPairs] = None,
    output_pairs_file: Optional[Path] = None,
    detailed: bool = False,
):
    """
    Write the gene with TA annotation to a file.
    """

    cluster_objects: Dict[int, GeneCluster] = {}  # cluster_id → GeneCluster
    single_genes_count = 0
    gene_in_cluster_count = 0
    cluster_score_cache: Dict[int, Dict[str, Any]] = {}  # cluster_id → score dict
    _score_null = (
        _SCORE_NULL_FULL
        if detailed
        else {k: v for k, v in _SCORE_NULL_FULL.items() if k not in _SCORE_DETAIL_KEYS}
    )
    with open(output_file, "w") as fl:
        writer = None
        for contig_name, cdss in contig_name_and_cdss_with_ta_hit:
            for cds in cdss:

                ta_gene_info: Dict[str, Any] = {}

                ta_gene_info["contig_name"] = contig_name
                ta_gene_info["gene_id"] = cds.id
                ta_gene_info["start"] = cds.start
                ta_gene_info["end"] = cds.stop
                ta_gene_info["strand"] = cds.strand.value
                ta_gene_info["length_aa"] = (cds.stop - cds.start + 1) // 3
                ta_gene_info["product"] = cds.product

                if cds.ta_cluster:
                    ta_gene_info["ta_system_id"] = cds.ta_cluster.id
                    ta_gene_info["is_single_gene"] = False
                    cluster_objects[cds.ta_cluster.id] = cds.ta_cluster
                    gene_in_cluster_count += 1
                else:
                    ta_gene_info["is_single_gene"] = True
                    ta_gene_info["ta_system_id"] = None
                    single_genes_count += 1

                hit_info = (
                    summarise_ta_hits(cds.ta_hits, hmm_db_info)
                    if detailed
                    else _summarise_ta_hits_simple(cds.ta_hits, hmm_db_info)
                )

                ta_gene_info.update(hit_info)

                if ref_stats is not None:
                    if cds.ta_cluster:
                        cid = cds.ta_cluster.id
                        if cid not in cluster_score_cache:
                            cluster_score_cache[cid] = score_cluster(
                                cds.ta_cluster,
                                hmm_db_info,
                                ref_stats,
                                known_pairs=known_pairs,
                            )
                        score_dict = cluster_score_cache[cid]
                        if not detailed:
                            score_dict = {
                                k: v
                                for k, v in score_dict.items()
                                if k not in _SCORE_DETAIL_KEYS
                            }
                        ta_gene_info.update(score_dict)
                    else:
                        ta_gene_info.update(_score_null)

                if writer is None:
                    writer = csv.DictWriter(
                        fl, fieldnames=ta_gene_info.keys(), delimiter="\t"
                    )
                    writer.writeheader()

                writer.writerow(ta_gene_info)

    logger.info(f"Finished writing genes with TA annotations to file '{output_file}'")
    logger.info(
        f"Identified {len(cluster_objects)} gene groups with TA annotations, "
        f"containing a total of {gene_in_cluster_count} genes."
    )
    logger.info(f"Identified {single_genes_count} single genes with TA annotations.")

    if output_pairs_file is not None:
        write_pairs_with_ta_annotation(
            cluster_objects,
            hmm_db_info,
            output_pairs_file,
            ref_stats=ref_stats,
            known_pairs=known_pairs,
            detailed=detailed,
        )


def write_pairs_with_ta_annotation(
    clusters: Dict[int, GeneCluster],
    hmm_db_info: Dict[str, Dict[str, str]],
    output_file: Path,
    ref_stats: Optional[RefStats] = None,
    known_pairs: Optional[KnownPairs] = None,
    detailed: bool = False,
):
    """Write one TSV row per valid (toxin, antitoxin) pair across all clusters.

    For clusters with more than one toxin or antitoxin, all combinations are
    written as individual rows, each independently scored.
    Clusters with no identifiable toxin or antitoxin gene are skipped.
    """
    pairs_written = 0
    skipped = 0

    with open(output_file, "w") as fl:
        writer = None

        for cluster_id, cluster in sorted(clusters.items()):
            genes = list(cluster.genes)
            types = [get_gene_type(g, hmm_db_info) for g in genes]

            toxins = [g for g, t in zip(genes, types) if t == "Toxin"]
            antitoxins = [g for g, t in zip(genes, types) if t == "Antitoxin"]

            if not toxins or not antitoxins:
                skipped += 1
                continue

            for toxin, antitoxin in itertools.product(toxins, antitoxins):
                # Contig (both genes share the same contig in a cluster)
                contig_name = toxin.contig_id

                # Structural features
                toxin_length_aa = (toxin.stop - toxin.start + 1) // 3
                antitoxin_length_aa = (antitoxin.stop - antitoxin.start + 1) // 3

                upstream = toxin if toxin.start <= antitoxin.start else antitoxin
                downstream = antitoxin if toxin.start <= antitoxin.start else toxin
                intergenic_distance = downstream.start - upstream.stop

                # HMM hit summaries per gene
                _hit_fn = summarise_ta_hits if detailed else _summarise_ta_hits_simple
                toxin_hits = _hit_fn(toxin.ta_hits, hmm_db_info)
                antitoxin_hits = _hit_fn(antitoxin.ta_hits, hmm_db_info)

                row: Dict[str, Any] = {
                    "ta_system_id": cluster_id,
                    "contig_name": contig_name,
                    "toxin_gene_id": toxin.id,
                    "toxin_start": toxin.start,
                    "toxin_end": toxin.stop,
                    "toxin_strand": toxin.strand.value,
                    "toxin_product": toxin.product,
                    "toxin_length_aa": toxin_length_aa,
                }
                row.update({f"toxin_{k}": v for k, v in toxin_hits.items()})

                row["antitoxin_gene_id"] = antitoxin.id
                row["antitoxin_start"] = antitoxin.start
                row["antitoxin_end"] = antitoxin.stop
                row["antitoxin_strand"] = antitoxin.strand.value
                row["antitoxin_product"] = antitoxin.product
                row["antitoxin_length_aa"] = antitoxin_length_aa
                row.update({f"antitoxin_{k}": v for k, v in antitoxin_hits.items()})

                row["intergenic_distance"] = intergenic_distance

                if ref_stats is not None:
                    score_info = score_pair(
                        toxin,
                        antitoxin,
                        hmm_db_info,
                        ref_stats,
                        known_pairs=known_pairs,
                    )
                    if not detailed:
                        score_info = {
                            k: v
                            for k, v in score_info.items()
                            if k not in _SCORE_DETAIL_KEYS
                        }
                    row.update(score_info)
                else:
                    _pair_score_null = (
                        _SCORE_NULL_FULL
                        if detailed
                        else {
                            k: v
                            for k, v in _SCORE_NULL_FULL.items()
                            if k not in _SCORE_DETAIL_KEYS
                        }
                    )
                    row.update(_pair_score_null)

                if not detailed:
                    row = {k: v for k, v in row.items() if k not in _PAIRS_DETAIL_KEYS}

                if writer is None:
                    writer = csv.DictWriter(fl, fieldnames=row.keys(), delimiter="\t")
                    writer.writeheader()

                writer.writerow(row)
                pairs_written += 1

    logger.info(
        f"Finished writing TA pairs to '{output_file}' ({pairs_written} pairs written)."
    )
    if skipped:
        logger.info(
            f"Skipped {skipped} clusters with no identifiable toxin or antitoxin gene."
        )
