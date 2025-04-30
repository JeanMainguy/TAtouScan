from typing import List, Generator, Any
import logging

from tatouscan.models import Cds


logger = logging.getLogger(__name__)


def group_cdss_with_ta_annotation(
    contig_name_and_cdss: Generator[tuple[str, List[Cds]], Any, None], max_distance: int
):
    for contig_name, cdss in contig_name_and_cdss:

        sorted_cdss_with_ta_hit = sorted([cds for cds in cdss if cds.ta_hits], key=lambda x: x.start)

        for i, cds_i in enumerate(sorted_cdss_with_ta_hit):

            for cds_j in sorted_cdss_with_ta_hit[i+1:]:

                if cds_i.distance_from(cds_j) <= max_distance:

                    cds_i.add_neigbor_gene(cds_j)

                    cds_i.update_ta_cluster(cds_j)

                else:
                    break

        yield contig_name, sorted_cdss_with_ta_hit
