from typing import Dict

from tatouscan.models import TaHit
from tatouscan.writer import (
    PFAM_INFO_COLUMNS,
    _summarise_ta_hits_simple,
    get_best_hit,
    get_pfam_info,
    summarise_ta_hits,
)

HMM_DB_INFO_WITH_PFAM: Dict[str, Dict[str, str]] = {
    "T.51159": {
        "type": "Toxin",
        "source": "TASmania",
        "supplementary_info": "YoeB-like toxin",
        "pfam_accession": "PF05016.22",
        "pfam_name": "ParE_toxin",
        "pfam_description": "ParE toxin of type II toxin-antitoxin system, parDE",
        "pfam_evalue": "9.4e-10",
        "pfam_coverage": "0.95",
        "pfam_hits": "ParE_toxin(PF05016.22):9.4e-10:5-87;YoeB_toxin(PF06769.21):6.8e-09:3-86",
    },
    "T.999": {  # profile without significant Pfam hit
        "type": "Toxin",
        "source": "TADB3",
        "supplementary_info": "None",
        "pfam_accession": "",
        "pfam_name": "",
        "pfam_description": "",
        "pfam_evalue": "",
        "pfam_coverage": "",
        "pfam_hits": "",
    },
}


def test_get_pfam_info_with_annotation():
    pfam = get_pfam_info(HMM_DB_INFO_WITH_PFAM["T.51159"])
    assert pfam == {
        "hmm_pfam_accession": "PF05016.22",
        "hmm_pfam_name": "ParE_toxin",
        "hmm_pfam_description": "ParE toxin of type II toxin-antitoxin system, parDE",
        "hmm_pfam_hits": "ParE_toxin(PF05016.22):9.4e-10:5-87;YoeB_toxin(PF06769.21):6.8e-09:3-86",
    }


def test_get_pfam_info_empty_values_are_none():
    assert get_pfam_info(HMM_DB_INFO_WITH_PFAM["T.999"]) == {
        col: None for col in PFAM_INFO_COLUMNS
    }


def test_summarise_simple_includes_pfam_columns():
    hits = [
        TaHit("prot1", "T.51159", score=50.0, evalue=1e-12),
        TaHit("prot1", "T.999", score=20.0, evalue=1e-3),
    ]
    summary = _summarise_ta_hits_simple(hits, HMM_DB_INFO_WITH_PFAM)
    assert summary["hmm_name"] == "T.51159"
    assert summary["gene_type"] == "Toxin"
    assert summary["hmm_pfam_name"] == "ParE_toxin"
    assert summary["hmm_pfam_accession"] == "PF05016.22"
    assert summary["hmm_pfam_hits"].startswith("ParE_toxin(PF05016.22)")


def test_get_best_hit_includes_pfam_columns():
    best = get_best_hit([TaHit("prot1", "T.51159", 50.0, 1e-12)], HMM_DB_INFO_WITH_PFAM)
    assert best["hmm_description"] == "YoeB-like toxin"
    assert best["hmm_pfam_name"] == "ParE_toxin"

    empty = get_best_hit([], HMM_DB_INFO_WITH_PFAM)
    assert set(PFAM_INFO_COLUMNS) <= set(empty)
    assert all(v is None for v in empty.values())


def test_summarise_detailed_profile_missing_from_hmm_info():
    """A profile present in ta.hmm but absent from hmm_info.tsv must not crash."""
    hits = [TaHit("prot1", "T.unknown", score=50.0, evalue=1e-12)]
    summary = summarise_ta_hits(hits, HMM_DB_INFO_WITH_PFAM)
    assert summary["gene_type"] is None
    assert "TASmania_hmm_pfam_name" in summary
