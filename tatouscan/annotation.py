from typing import List
import pyhmmer.easel

from pyhmmer.plan7 import HMM
from collections import defaultdict

from pathlib import Path

from tatouscan.models import TaHit


def find_ta_hits(faa_file: Path, hmm_db: Path, e_value_threshold: float = 0.01):
    """
    Find hits of TAs in a protein fasta file.

    :params faa_file: Path to a protein fasta file.
    :params hmm_db: Path to a HMM database.
    :params e_value_threshold: E-value threshold for hits. Default is 0.01.

    :returns: A dictionary with protein ids as keys and a list of TaHit objects as values.
    """
    proteins = None
    with pyhmmer.easel.SequenceFile(faa_file, digital=True) as seqs_file:
        proteins = seqs_file.read_block()

    if proteins is None:
        raise ValueError(f"No sequences found in {faa_file}")

    hmms: List[HMM] = []
    with pyhmmer.plan7.HMMFile(hmm_db) as hmm_file:
        for hmm in hmm_file:
            hmms.append(hmm)

    protein_id_to_hits: defaultdict[str, List[TaHit]] = defaultdict(list)

    for hits in pyhmmer.hmmsearch(hmms, proteins, E=e_value_threshold):
        ta_name = hits.query.name.decode()

        for hit in hits:
            protein_id = hit.name.decode()
            protein_id_to_hits[protein_id].append(
                TaHit(protein_id, ta_name, hit.score, hit.evalue)
            )

    return protein_id_to_hits
