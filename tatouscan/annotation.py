from typing import List
import pyhmmer.easel
import logging
from rich.progress import track


from pyhmmer.plan7 import HMM
from collections import defaultdict

from pathlib import Path

from tatouscan.models import TaHit

logger = logging.getLogger(__name__)


def find_ta_hits(faa_file: Path, hmm_db: Path, e_value_threshold: float = 0.01):
    """
    Find hits of TAs in a protein fasta file.

    :params faa_file: Path to a protein fasta file.
    :params hmm_db: Path to a HMM database.
    :params e_value_threshold: E-value threshold for hits. Default is 0.01.

    :returns: A dictionary with protein ids as keys and a list of TaHit objects as values.
    """
    logger.info(f"Loading HMM profiles from database: '{hmm_db}'")
    hmms: List[HMM] = []
    with pyhmmer.plan7.HMMFile(hmm_db) as hmm_file:
        for hmm in hmm_file:
            hmms.append(hmm)
    logger.info(f"Loaded {len(hmms)} HMM profiles from '{hmm_db}'")

    logger.info(f"Reading protein sequences from FASTA file: '{faa_file}'")
    with pyhmmer.easel.SequenceFile(faa_file, digital=True) as seqs_file:
        proteins = seqs_file.read_block()

        logger.info(f"Found {len(proteins)} protein sequences in '{faa_file}'")

        logger.info(
            f"Starting TA hit search with E-value threshold: {e_value_threshold}"
        )

        protein_id_to_hits: defaultdict[str, List[TaHit]] = defaultdict(list)

        for hits in track(
            pyhmmer.hmmsearch(hmms, proteins, E=e_value_threshold),
            total=len(proteins),
            description="Searching for TA hits",
        ):
            ta_name = hits.query.name.decode()

            for hit in hits:
                protein_id = hit.name.decode()

                protein_id_to_hits[protein_id].append(
                    TaHit(protein_id, ta_name, hit.score, hit.evalue)
                )

        logger.info(f"Identified {len(protein_id_to_hits)} proteins with TA hits")
        logger.info("TA hit search completed successfully")

        return protein_id_to_hits
