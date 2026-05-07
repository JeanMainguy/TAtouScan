"""
database.py – TAtouScan reference database loader.

A TAtouScan database is a directory containing four fixed files:

    ta.hmm                  HMMER3 HMM profiles (toxin + antitoxin)
    hmm_info.tsv            HMM metadata (name, type, source, …)
    family_statistics.tsv   Per-family reference statistics for scoring
    known_pairs.tsv         Known (toxin_family, AT_family) co-occurrence pairs

The directory is typically distributed as a .tar.gz archive that users
extract once before running TAtouScan.  Pass the extracted directory to
the ``--db`` CLI option.

Example::

    tatouscan --db /path/to/tatouscan_db/ --gff genome.gff --faa genome.faa …
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

# Fixed filenames inside the database directory.
DB_FILENAMES: dict[str, str] = {
    "hmm_db": "ta.hmm",
    "hmm_info": "hmm_info.tsv",
    "ref_stats": "family_statistics.tsv",
    "known_pairs": "known_pairs.tsv",
}


@dataclass(frozen=True)
class TAtouScanDB:
    """Resolved paths to the four components of a TAtouScan database directory."""

    hmm_db: Path
    hmm_info: Path
    ref_stats: Path
    known_pairs: Path


def load_db(path: Path) -> TAtouScanDB:
    """Validate *path* as a TAtouScan database directory and return resolved paths.

    Raises
    ------
    NotADirectoryError
        If *path* does not exist or is not a directory.
    FileNotFoundError
        If any of the four required files is missing from the directory.
    """
    if not path.exists():
        raise NotADirectoryError(f"Database path does not exist: {path}")
    if not path.is_dir():
        raise NotADirectoryError(
            f"--db must point to a directory (not a file): {path}\n"
            "If you have a .tar.gz archive, extract it first: tar -xzf <archive>"
        )

    missing = [fname for fname in DB_FILENAMES.values() if not (path / fname).exists()]
    if missing:
        raise FileNotFoundError(
            f"TAtouScan database at '{path}' is missing required file(s): "
            + ", ".join(missing)
            + f"\nExpected: {', '.join(DB_FILENAMES.values())}"
        )

    return TAtouScanDB(
        hmm_db=path / DB_FILENAMES["hmm_db"],
        hmm_info=path / DB_FILENAMES["hmm_info"],
        ref_stats=path / DB_FILENAMES["ref_stats"],
        known_pairs=path / DB_FILENAMES["known_pairs"],
    )
