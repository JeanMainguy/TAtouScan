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


# Pfam annotation of the HMM profiles (output columns, in output order).
# Empty values (profile without significant Pfam hit) are written as None.
#   hmm_pfam_name         closest Pfam family            (hmm_info.tsv: pfam_name)
#   hmm_pfam_families     all significant family names   (derived from pfam_hits)
#   hmm_pfam_accession    closest family accession       (hmm_info.tsv: pfam_accession)
#   hmm_pfam_description  closest family description     (hmm_info.tsv: pfam_description)
#   hmm_pfam_hits         all significant families with  (hmm_info.tsv: pfam_hits)
#                         E-values and positions
PFAM_INFO_COLUMNS: tuple[str, ...] = (
    "hmm_pfam_name",
    "hmm_pfam_families",
    "hmm_pfam_accession",
    "hmm_pfam_description",
    "hmm_pfam_hits",
)

# Pfam columns written for the best hit in default output; the others are --detailed only.
PFAM_DEFAULT_COLUMNS: frozenset[str] = frozenset({"hmm_pfam_name", "hmm_pfam_families"})

# Pfam columns written for each per-source best hit (--detailed only).
PFAM_PER_SOURCE_COLUMNS: frozenset[str] = frozenset({"hmm_pfam_name"})


def pfam_family_names(pfam_hits: str) -> str:
    """Reduce ``NAME(ACC):E-value:positions;...`` to ``NAME;...``."""
    return ";".join(hit.split("(", 1)[0] for hit in pfam_hits.split(";") if hit)


def get_pfam_info(
    info: Dict[str, str], columns: Optional[frozenset[str]] = None
) -> Dict[str, Optional[str]]:
    """Return the Pfam annotation of an HMM profile from its hmm_info.tsv row.

    :param columns: output columns to keep (all Pfam columns if None).
    """
    pfam_hits = info.get("pfam_hits") or ""
    values = {
        "hmm_pfam_name": info.get("pfam_name"),
        "hmm_pfam_families": pfam_family_names(pfam_hits),
        "hmm_pfam_accession": info.get("pfam_accession"),
        "hmm_pfam_description": info.get("pfam_description"),
        "hmm_pfam_hits": pfam_hits,
    }
    return {
        col: values[col] or None
        for col in PFAM_INFO_COLUMNS
        if columns is None or col in columns
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
            **get_pfam_info({}, PFAM_PER_SOURCE_COLUMNS),
        }

    # Sort hits by e-value and return the first one
    best_hit = min(ta_hits, key=lambda hit: hit.evalue)
    info = hmm_db_info.get(best_hit.ta_name, {})

    return {
        "hmm_name": best_hit.ta_name,
        "hmm_score": best_hit.score,
        "hmm_evalue": best_hit.evalue,
        "hmm_description": info.get("supplementary_info"),
        **get_pfam_info(info, PFAM_PER_SOURCE_COLUMNS),
    }


def _summarise_ta_hits_simple(
    ta_hits: List[TaHit],
    hmm_db_info: Dict[str, Dict[str, str]],
    detailed: bool = False,
) -> Dict[str, Any]:
    """Return gene type and the single best HMM hit across all database sources.

    In detailed mode, all Pfam columns of the best hit are included; otherwise
    only the closest Pfam family name and the list of all significant Pfam hits.
    """
    best_hit = min(ta_hits, key=lambda hit: hit.evalue)
    info = hmm_db_info.get(best_hit.ta_name, {})
    return {
        "gene_type": info.get("type"),
        "hmm_name": best_hit.ta_name,
        "hmm_score": best_hit.score,
        "hmm_evalue": best_hit.evalue,
        "hmm_source": info.get("source"),
        "hmm_description": info.get("supplementary_info"),
        **get_pfam_info(info, None if detailed else PFAM_DEFAULT_COLUMNS),
    }


def summarise_ta_hits(ta_hits: List[TaHit], hmm_db_info: Dict[str, Dict[str, str]]):
    """Detailed hit summary: the best hit across all sources (with all Pfam
    columns), followed by the best hit of each source (TASmania, Other, TADB3)."""

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

    hit_info: Dict[str, Any] = _summarise_ta_hits_simple(
        ta_hits, hmm_db_info, detailed=True
    )

    hit_info.update({f"TASmania_{k}": v for k, v in best_tasmania_hit.items()})
    hit_info.update({f"Other_{k}": v for k, v in best_other_hit.items()})
    hit_info.update({f"TADB3_{k}": v for k, v in best_tadb3_hit.items()})

    return hit_info


def _hit_columns(detailed: bool) -> List[str]:
    """Columns of the HMM hit summary of a gene, in output order."""
    summarise = summarise_ta_hits if detailed else _summarise_ta_hits_simple
    return list(summarise([TaHit("", "", 0.0, 0.0)], {}).keys())


def _score_columns(detailed: bool) -> List[str]:
    """Columns of the score of a TA system or pair, in output order."""
    return [k for k in _SCORE_NULL_FULL if detailed or k not in _SCORE_DETAIL_KEYS]


def gene_output_columns(detailed: bool = False, with_score: bool = True) -> List[str]:
    """Columns of the gene output file, in output order."""
    columns = [
        "contig_name",
        "gene_id",
        "start",
        "end",
        "strand",
        "length_aa",
        "product",
        "ta_system_id",
        "is_single_gene",
    ]
    columns += _hit_columns(detailed)
    if with_score:
        columns += _score_columns(detailed)
    return columns


def pair_output_columns(detailed: bool = False) -> List[str]:
    """Columns of the pairs output file, in output order."""
    hit_columns = _hit_columns(detailed)
    columns = ["ta_system_id", "contig_name"]
    for role in ("toxin", "antitoxin"):
        columns += [
            f"{role}_{col}"
            for col in ("gene_id", "start", "end", "strand", "product", "length_aa")
        ]
        columns += [f"{role}_{col}" for col in hit_columns]
    columns.append("intergenic_distance")
    columns += _score_columns(detailed)
    if not detailed:
        columns = [col for col in columns if col not in _PAIRS_DETAIL_KEYS]
    return columns


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
    # The header is written even when no gene has a TA hit, so that the
    # output file can always be read as a table.
    with open(output_file, "w", newline="") as fl:
        writer = csv.DictWriter(
            fl,
            fieldnames=gene_output_columns(detailed, with_score=ref_stats is not None),
            delimiter="\t",
        )
        writer.writeheader()
        for contig_name, cdss in contig_name_and_cdss_with_ta_hit:
            for cds in cdss:

                ta_gene_info: Dict[str, Any] = {}

                ta_gene_info["contig_name"] = contig_name
                ta_gene_info["gene_id"] = cds.id
                ta_gene_info["start"] = cds.start
                ta_gene_info["end"] = cds.stop
                ta_gene_info["strand"] = cds.strand.value
                ta_gene_info["length_aa"] = cds.length_aa
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

    with open(output_file, "w", newline="") as fl:
        writer = csv.DictWriter(
            fl, fieldnames=pair_output_columns(detailed), delimiter="\t"
        )
        writer.writeheader()

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
                toxin_length_aa = toxin.length_aa
                antitoxin_length_aa = antitoxin.length_aa

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

                writer.writerow(row)
                pairs_written += 1

    logger.info(
        f"Finished writing TA pairs to '{output_file}' ({pairs_written} pairs written)."
    )
    if skipped:
        logger.info(
            f"Skipped {skipped} clusters with no identifiable toxin or antitoxin gene."
        )
