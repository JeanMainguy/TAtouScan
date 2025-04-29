from typing import List, Generator, Any
import logging
from rich.progress import track


from collections import defaultdict

from pathlib import Path

from tatouscan.models import TaHit, Cds


logger = logging.getLogger(__name__)


def group_cdss_with_ta_annotation(
    contig_name_and_cdss: Generator[tuple[str, List[Cds]], Any, None], max_distance: int
):
    print("Grouping CDSs with TA annotation...")
    for contig_name, cdss in contig_name_and_cdss:
        
        # ta_clusters:list[List[Cds]] = []
        sorted_cdss_with_ta_hit = sorted([cds for cds in cdss if cds.ta_hits], key=lambda x: x.start)

        for i, cds_i in enumerate(sorted_cdss_with_ta_hit):
            
            ta_cluster = [ cds_i ]

            for cds_j in sorted_cdss_with_ta_hit[i+1:]:
                
                if cds_i.distance_from(cds_j) <= max_distance:
                    
                    cds_i.add_neigbor_gene(cds_j)

                    cds_i.update_ta_cluster(cds_j)

                    ta_cluster.append(cds_j)

                else:
                    break
                                

        for cds in sorted_cdss_with_ta_hit:
            print(cds, cds.neighbor_genes)

        ta_clusters = {cds.ta_cluster for cds in sorted_cdss_with_ta_hit if cds.ta_cluster is not None}

        if ta_clusters:
            print(ta_clusters)

        input()