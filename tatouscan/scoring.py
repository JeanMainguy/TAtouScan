"""
scoring.py – Family-aware robust Z-score scoring of putative TA pairs.

For each two-gene cluster, compares three structural features
(toxin length, antitoxin length, intergenic distance) against the
distribution observed in known TADB3 systems.  Features are normalised
by family-specific reference statistics when available, falling back to
global statistics for rare or unrecognised families.

Robust Z-score (median + MAD)
------------------------------
    z_robust = (x − median) / (MAD / 0.6745)

The 0.6745 factor (= Φ⁻¹(0.75)) makes MAD equivalent to σ for normally
distributed data, so scores are directly comparable to classic Z-scores.
Unlike mean + std, median and MAD are unaffected by outliers and skewed
tails — which is important here because many TA families have non-normal
size distributions.

Unified score
-------------
    score = exp(−mean(|z_i|))

The mean is taken over all available z-score terms:
  • toxin_size_z, at_size_z, intergenic_distance_z  (structural features)
  • z_compat: 0 if the (toxin, AT) pair is known in TADB3,
              COMPAT_PENALTY (default 2.0) if the pair is unknown,
              excluded from the mean if family identity is unavailable.

Range: (0, 1].  1 = perfect match to the family median, known pair.
The exponential transform preserves ranking while giving an intuitive
[0, 1] scale for users.
"""

from __future__ import annotations

import csv
import itertools
import math
from pathlib import Path
from typing import Any, Dict, FrozenSet, List, Optional, Set, Tuple

from tatouscan.models import Cds, GeneCluster

# Features that are scored; names must match columns in the reference TSV
# and the keys returned in the score dict.
SCORED_FEATURES = ("toxin_size", "at_size", "intergenic_distance")

# Sentinel used for the global fallback row in the reference stats file
GLOBAL_FAMILY = "__global__"

# Family statistics type:
#   family → "n_pairs" → int
#   family → feature   → (center, scale)
#                         center = median (or mean for old TSVs)
#                         scale  = MAD/0.6745 (or std for old TSVs / MAD=0 fallback)
RefStats = Dict[str, Dict[str, Any]]

# HMM database metadata: hmm_name → {"type", "source", "supplementary_info", ...}
HmmDbInfo = Dict[str, Dict[str, str]]

# Known toxin–antitoxin family co-occurrence pairs from TADB3.
# Each element is (toxin_family, at_family).
KnownPairs = FrozenSet[Tuple[str, str]]

# Z-score equivalent penalty applied when a (toxin, AT) family pair is
# not found in the TADB3 known-pairs set.  2.0 means "an unknown pair is
# treated as a 2-sigma structural mismatch" in the unified score.
COMPAT_PENALTY = 2.0


# ---------------------------------------------------------------------------
# Loading reference data
# ---------------------------------------------------------------------------


def load_reference_statistics(path: Path) -> RefStats:
    """Load per-family reference statistics from *tadb3_family_statistics.tsv*.

    Prefers robust statistics (median + MAD) when the columns
    ``{feat}_median`` and ``{feat}_mad`` are present (produced by
    ``build_reference_features.py`` ≥ v2).  Falls back to mean + std for
    older TSV files that lack those columns.

    When MAD = 0 (e.g. a family where all representatives have identical
    values), falls back to std for that feature so that ``_z_score`` can
    still return a useful value rather than *None*.

    Returns a dict:  family → feature → (center, scale)
    Where center = median (or mean) and scale = MAD/0.6745 (or std).
    The special key ``"__global__"`` holds the global-fallback stats.
    """
    MAD_CONSISTENCY = 0.6745  # Φ⁻¹(0.75): makes MAD ≡ σ for normal data

    stats: RefStats = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames or []
        has_robust = all(
            f"{feat}_median" in fieldnames and f"{feat}_mad" in fieldnames
            for feat in SCORED_FEATURES
        )
        for row in reader:
            family = row["family"]
            entry: Dict[str, Any] = {}
            for feat in SCORED_FEATURES:
                if has_robust:
                    center = float(row[f"{feat}_median"])
                    mad = float(row[f"{feat}_mad"])
                    scale = mad / MAD_CONSISTENCY
                    if scale == 0.0:
                        # MAD=0: all values identical — fall back to std
                        scale = float(row.get(f"{feat}_std", 0.0))
                else:
                    center = float(row[f"{feat}_mean"])
                    scale = float(row[f"{feat}_std"])
                entry[feat] = (center, scale)
            entry["n_pairs"] = int(row["n_pairs"])
            stats[family] = entry
    return stats


def load_known_pairs(path: Path) -> KnownPairs:
    """Load the set of known (toxin_family, at_family) co-occurrence pairs.

    Reads *tadb3_known_pairs.tsv* produced by build_compatibility_matrix.py.
    Each row must have at minimum the columns ``toxin_family`` and ``at_family``.
    """
    pairs: Set[Tuple[str, str]] = set()
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pairs.add((row["toxin_family"], row["at_family"]))
    return frozenset(pairs)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def get_gene_type(cds: Cds, hmm_db_info: HmmDbInfo) -> Optional[str]:
    """Return the type string (``"Toxin"`` or ``"Antitoxin"``) for *cds*."""
    if not cds.ta_hits:
        return None
    best = min(cds.ta_hits, key=lambda h: h.evalue)
    info = hmm_db_info.get(best.ta_name)
    if info is None:
        return None
    return info.get("type")  # "Toxin" or "Antitoxin"


def _best_tadb3_toxin_family(cds: Cds, hmm_db_info: HmmDbInfo) -> Optional[str]:
    """Return the TADB3 toxin HMM profile name (= family ID) for *cds*.

    Only considers hits against TADB3 toxin profiles.  Returns *None* if
    the gene has no TADB3 toxin hit.
    """
    tadb3_toxin_hits = [
        h
        for h in cds.ta_hits
        if h.ta_name in hmm_db_info
        and hmm_db_info[h.ta_name].get("source") == "TADB3"
        and hmm_db_info[h.ta_name].get("type") == "Toxin"
    ]
    if not tadb3_toxin_hits:
        return None
    return min(tadb3_toxin_hits, key=lambda h: h.evalue).ta_name


def _best_tadb3_at_family(cds: Cds, hmm_db_info: HmmDbInfo) -> Optional[str]:
    """Return the TADB3 antitoxin HMM profile name (= family ID) for *cds*.

    Only considers hits against TADB3 antitoxin profiles.  Returns *None* if
    the gene has no TADB3 antitoxin hit.
    """
    tadb3_at_hits = [
        h
        for h in cds.ta_hits
        if h.ta_name in hmm_db_info
        and hmm_db_info[h.ta_name].get("source") == "TADB3"
        and hmm_db_info[h.ta_name].get("type") == "Antitoxin"
    ]
    if not tadb3_at_hits:
        return None
    return min(tadb3_at_hits, key=lambda h: h.evalue).ta_name


def _z_score(value: float, mean: float, std: float) -> Optional[float]:
    if std == 0.0:
        return None
    return (value - mean) / std


def _round_optional(value: Optional[float], ndigits: int = 3) -> Optional[float]:
    return None if value is None else round(value, ndigits)


# ---------------------------------------------------------------------------
# Cluster scoring
# ---------------------------------------------------------------------------


def score_cluster(
    cluster: GeneCluster,
    hmm_db_info: HmmDbInfo,
    ref_stats: RefStats,
    known_pairs: Optional[KnownPairs] = None,
) -> Dict[str, Any]:
    """Compute family-aware robust Z-score features for a two-gene *cluster*.

    Returns a dict with keys:
        matched_family          – family used for look-up (or ``"__global__"``)
        n_reference_pairs       – number of reference pairs in that family
        toxin_size_z            – robust Z-score of toxin amino-acid length
        at_size_z               – robust Z-score of antitoxin amino-acid length
        intergenic_distance_z   – robust Z-score of intergenic distance
        pair_is_known           – 1/0 if (toxin, AT) pair is in TADB3; None if
                                  family identity is unavailable or known_pairs
                                  was not provided
        score                   – unified match score in (0, 1]:
                                  ``exp(−mean(|z_structural| + z_compat))``
                                  where z_compat = 0 (known), COMPAT_PENALTY
                                  (unknown), or excluded (None)

    Scores use median + MAD internally: ``z = (x − median) / (MAD / 0.6745)``.
    This is robust to the non-normal (right-skewed) size distributions that
    are common across TA families.  If MAD = 0 for a feature, std is used
    instead.  If std = 0 as well, that feature returns *None*.

    All numeric fields are *None* when scoring cannot proceed (cluster has
    no identifiable toxin+antitoxin pair).  For clusters with more than one
    toxin or antitoxin, all (toxin, antitoxin) combinations are scored and
    the best-scoring pair is returned.
    """
    null_result: Dict[str, Any] = {
        "matched_family": None,
        "n_reference_pairs": None,
        "toxin_size_z": None,
        "at_size_z": None,
        "intergenic_distance_z": None,
        "pair_is_known": None,
        "score": None,
    }

    genes = list(cluster.genes)
    types = [get_gene_type(g, hmm_db_info) for g in genes]
    toxins: List[Cds] = [g for g, t in zip(genes, types) if t == "Toxin"]
    antitoxins: List[Cds] = [g for g, t in zip(genes, types) if t == "Antitoxin"]

    if not toxins or not antitoxins:
        return null_result

    best: Optional[Dict[str, Any]] = None
    for toxin, antitoxin in itertools.product(toxins, antitoxins):
        result = score_pair(toxin, antitoxin, hmm_db_info, ref_stats, known_pairs)
        if best is None or (
            result["score"] is not None
            and (best["score"] is None or result["score"] > best["score"])
        ):
            best = result
    return best if best is not None else null_result


def score_pair(
    toxin: Cds,
    antitoxin: Cds,
    hmm_db_info: HmmDbInfo,
    ref_stats: RefStats,
    known_pairs: Optional[KnownPairs] = None,
) -> Dict[str, Any]:
    """Score a single (toxin, antitoxin) gene pair."""
    # Protein length in amino acids from 1-based inclusive nucleotide coordinates
    toxin_size = (toxin.stop - toxin.start + 1) // 3
    at_size = (antitoxin.stop - antitoxin.start + 1) // 3

    # Intergenic distance (negative = overlap), matching reference formula:
    # downstream.genomic_left − upstream.genomic_right
    if toxin.start <= antitoxin.start:
        upstream, downstream = toxin, antitoxin
    else:
        upstream, downstream = antitoxin, toxin
    intergenic_distance = downstream.start - upstream.stop

    # Family look-up: use toxin's best TADB3 hit as family ID
    family = _best_tadb3_toxin_family(toxin, hmm_db_info)
    if not family or family not in ref_stats:
        family = GLOBAL_FAMILY

    fam_stats: Dict[str, Any] = ref_stats.get(family, {})

    # Z-scores
    def _z(value: float, feat: str) -> Optional[float]:
        if feat not in fam_stats:
            return None
        mean, std = fam_stats[feat]
        return _z_score(value, mean, std)

    z_toxin = _z(float(toxin_size), "toxin_size")
    z_at = _z(float(at_size), "at_size")
    z_intergenic = _z(float(intergenic_distance), "intergenic_distance")

    # Compatibility term: is this (toxin_family, at_family) a known TADB3 pair?
    pair_is_known: Optional[int] = None
    z_compat: Optional[float] = None
    if known_pairs is not None:
        toxin_family_id = _best_tadb3_toxin_family(toxin, hmm_db_info)
        at_family_id = _best_tadb3_at_family(antitoxin, hmm_db_info)
        if toxin_family_id is not None and at_family_id is not None:
            pair_is_known = 1 if (toxin_family_id, at_family_id) in known_pairs else 0
            z_compat = 0.0 if pair_is_known == 1 else COMPAT_PENALTY

    # Unified score: exp(−mean(|z_i|)) over all available terms
    # z_compat is excluded when family identity is unavailable (stays None)
    z_terms = [abs(z) for z in (z_toxin, z_at, z_intergenic) if z is not None]
    if z_compat is not None:
        z_terms.append(z_compat)
    score: Optional[float] = (
        round(math.exp(-sum(z_terms) / len(z_terms)), 4) if z_terms else None
    )

    return {
        "matched_family": family,
        "n_reference_pairs": fam_stats.get("n_pairs"),
        "toxin_size_z": _round_optional(z_toxin),
        "at_size_z": _round_optional(z_at),
        "intergenic_distance_z": _round_optional(z_intergenic),
        "pair_is_known": pair_is_known,
        "score": score,
    }
