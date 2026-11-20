from typing import Dict

import pytest

from tatouscan.models import TaHit
from tatouscan.writer import (
    PFAM_DEFAULT_COLUMNS,
    PFAM_INFO_COLUMNS,
    PFAM_PER_SOURCE_COLUMNS,
    _summarise_ta_hits_simple,
    gene_output_columns,
    get_best_hit,
    get_pfam_info,
    pair_output_columns,
    pfam_family_names,
    summarise_ta_hits,
    write_gene_with_ta_annotation,
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
        "hmm_pfam_name": "ParE_toxin",
        "hmm_pfam_families": "ParE_toxin;YoeB_toxin",
        "hmm_pfam_accession": "PF05016.22",
        "hmm_pfam_description": "ParE toxin of type II toxin-antitoxin system, parDE",
        "hmm_pfam_hits": "ParE_toxin(PF05016.22):9.4e-10:5-87;YoeB_toxin(PF06769.21):6.8e-09:3-86",
    }
    assert list(pfam) == list(PFAM_INFO_COLUMNS)


def test_pfam_family_names():
    assert pfam_family_names("") == ""
    assert pfam_family_names("PIN(PF01850.28):5.4e-24:101-220") == "PIN"
    assert (
        pfam_family_names(
            "TPR_12(PF13424.13):3.4e-15:649-716,732-797;HTH_5(PF01022.27):5.5e-19:266-311"
        )
        == "TPR_12;HTH_5"
    )


def test_get_pfam_info_empty_values_are_none():
    assert get_pfam_info(HMM_DB_INFO_WITH_PFAM["T.999"]) == {
        col: None for col in PFAM_INFO_COLUMNS
    }


def test_get_pfam_info_column_subset():
    assert get_pfam_info(HMM_DB_INFO_WITH_PFAM["T.51159"], PFAM_DEFAULT_COLUMNS) == {
        "hmm_pfam_name": "ParE_toxin",
        "hmm_pfam_families": "ParE_toxin;YoeB_toxin",
    }


HITS = [
    TaHit("prot1", "T.51159", score=50.0, evalue=1e-12),
    TaHit("prot1", "T.999", score=20.0, evalue=1e-3),
]


def test_summarise_simple_default_pfam_columns():
    """Default output: only the closest Pfam name and the names of all Pfam families."""
    summary = _summarise_ta_hits_simple(HITS, HMM_DB_INFO_WITH_PFAM)
    assert summary["hmm_name"] == "T.51159"
    assert summary["gene_type"] == "Toxin"
    assert summary["hmm_pfam_name"] == "ParE_toxin"
    assert summary["hmm_pfam_families"] == "ParE_toxin;YoeB_toxin"
    for detail_col in ("hmm_pfam_accession", "hmm_pfam_description", "hmm_pfam_hits"):
        assert detail_col not in summary


def test_summarise_simple_detailed_pfam_columns():
    summary = _summarise_ta_hits_simple(HITS, HMM_DB_INFO_WITH_PFAM, detailed=True)
    assert set(PFAM_INFO_COLUMNS) <= set(summary)
    assert summary["hmm_pfam_accession"] == "PF05016.22"
    assert summary["hmm_pfam_hits"].startswith("ParE_toxin(PF05016.22):9.4e-10")


def test_summarise_detailed_is_superset_of_default():
    """--detailed keeps every default column and adds per-source columns."""
    default = _summarise_ta_hits_simple(HITS, HMM_DB_INFO_WITH_PFAM)
    detailed = summarise_ta_hits(HITS, HMM_DB_INFO_WITH_PFAM)
    assert set(default) <= set(detailed)
    assert {k: detailed[k] for k in default} == default
    assert detailed["hmm_pfam_description"].startswith("ParE toxin")
    assert detailed["TASmania_hmm_name"] == "T.51159"
    assert detailed["TASmania_hmm_pfam_name"] == "ParE_toxin"
    assert detailed["TADB3_hmm_name"] == "T.999"
    assert detailed["TADB3_hmm_pfam_name"] is None


def test_get_best_hit_per_source_pfam_columns():
    """Per-source best hits carry only the Pfam family name."""
    best = get_best_hit([TaHit("prot1", "T.51159", 50.0, 1e-12)], HMM_DB_INFO_WITH_PFAM)
    assert best["hmm_description"] == "YoeB-like toxin"
    assert best["hmm_pfam_name"] == "ParE_toxin"
    assert not set(PFAM_INFO_COLUMNS) - PFAM_PER_SOURCE_COLUMNS & set(best)

    empty = get_best_hit([], HMM_DB_INFO_WITH_PFAM)
    assert set(empty) == set(best)
    assert all(v is None for v in empty.values())


def test_summarise_detailed_profile_missing_from_hmm_info():
    """A profile present in ta.hmm but absent from hmm_info.tsv must not crash."""
    hits = [TaHit("prot1", "T.unknown", score=50.0, evalue=1e-12)]
    summary = summarise_ta_hits(hits, HMM_DB_INFO_WITH_PFAM)
    assert summary["gene_type"] is None
    assert "TASmania_hmm_pfam_name" in summary


def _read_tsv(path):
    import csv

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return reader.fieldnames, list(reader)


def _toxin_antitoxin_pair():
    from tatouscan.models import Cds, Frame, Strand

    toxin = Cds("tox", "contig1", [(1, 300)], Strand("+"), Frame(0))
    toxin.ta_hits = [TaHit("tox", "T.51159", 100.0, 1e-20)]
    antitoxin = Cds("at", "contig1", [(301, 540)], Strand("+"), Frame(0))
    antitoxin.ta_hits = [TaHit("at", "AT.1", 80.0, 1e-15)]
    toxin.update_ta_cluster(antitoxin)
    return toxin, antitoxin


HMM_DB_INFO_PAIR: Dict[str, Dict[str, str]] = {
    **HMM_DB_INFO_WITH_PFAM,
    "AT.1": {"type": "Antitoxin", "source": "TASmania", "supplementary_info": ""},
}


@pytest.mark.parametrize("detailed", [False, True])
def test_write_outputs_without_hits_have_header(tmp_path, detailed):
    gene_file, pairs_file = tmp_path / "genes.tsv", tmp_path / "pairs.tsv"

    write_gene_with_ta_annotation(
        iter([("contig1", [])]),
        HMM_DB_INFO_PAIR,
        gene_file,
        ref_stats={},
        known_pairs=frozenset(),
        output_pairs_file=pairs_file,
        detailed=detailed,
    )

    assert _read_tsv(gene_file) == (gene_output_columns(detailed), [])
    assert _read_tsv(pairs_file) == (pair_output_columns(detailed), [])


@pytest.mark.parametrize("detailed", [False, True])
def test_write_outputs_header_matches_rows(tmp_path, detailed):
    gene_file, pairs_file = tmp_path / "genes.tsv", tmp_path / "pairs.tsv"
    toxin, antitoxin = _toxin_antitoxin_pair()

    write_gene_with_ta_annotation(
        iter([("contig1", [toxin, antitoxin])]),
        HMM_DB_INFO_PAIR,
        gene_file,
        ref_stats={},
        known_pairs=frozenset(),
        output_pairs_file=pairs_file,
        detailed=detailed,
    )

    gene_columns, gene_rows = _read_tsv(gene_file)
    assert gene_columns == gene_output_columns(detailed)
    assert [row["gene_id"] for row in gene_rows] == ["tox", "at"]
    assert all(row["ta_system_id"] and row["hmm_name"] for row in gene_rows)

    pair_columns, pair_rows = _read_tsv(pairs_file)
    assert pair_columns == pair_output_columns(detailed)
    assert len(pair_rows) == 1
    assert pair_rows[0]["toxin_hmm_name"] == "T.51159"
    assert pair_rows[0]["antitoxin_hmm_name"] == "AT.1"
    assert pair_rows[0]["intergenic_distance"] == "1"
