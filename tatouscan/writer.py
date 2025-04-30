from typing import Any, Generator, Dict, List

from tatouscan.models import Cds, TaHit
from pathlib import Path
import csv
import logging

logger = logging.getLogger(__name__)


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


def summarise_ta_hits(ta_hits: List[TaHit], hmm_db_info: Dict[str, Dict[str, str]]):

    # gene_types = [
    #     hmm_db_info[ta_hit.ta_name]["type"]
    #     for ta_hit in ta_hits
    #     if ta_hit.ta_name in hmm_db_info
    # ]
    # if gene_types and gene_types.count("Toxin") / len(gene_types) > 0.75:
    #     gene_type = "Toxin"
    # elif gene_types and gene_types.count("Antitoxin") / len(gene_types) > 0.75:
    #     gene_type = "Antitoxin"
    # else:
    #     gene_type = "Unknown"

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
):
    """
    Write the gene with TA annotation to a file.
    """
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
                ta_gene_info["product"] = cds.product

                if cds.ta_cluster:
                    ta_gene_info["ta_system_id"] = cds.ta_cluster.id
                    ta_gene_info["is_single_gene"] = False
                else:
                    ta_gene_info["is_single_gene"] = True
                    ta_gene_info["ta_system_id"] = None

                hit_info = summarise_ta_hits(cds.ta_hits, hmm_db_info)

                ta_gene_info.update(hit_info)
                logger.info(ta_gene_info)

                if writer is None:
                    print(ta_gene_info)
                    writer = csv.DictWriter(
                        fl, fieldnames=ta_gene_info.keys(), delimiter="\t"
                    )
                    writer.writeheader()

                writer.writerow(ta_gene_info)
