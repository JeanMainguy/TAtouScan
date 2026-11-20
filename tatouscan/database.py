"""
database.py – TAtouScan reference database loader.

A TAtouScan database is a directory containing five fixed files:

    version.txt             Database version (e.g. ``0.3.0``)
    ta.hmm                  HMMER3 HMM profiles (toxin + antitoxin)
    hmm_info.tsv            HMM metadata (name, type, source, Pfam annotation, …)
    family_statistics.tsv   Per-family reference statistics for scoring
    known_pairs.tsv         Known (toxin_family, AT_family) co-occurrence pairs

The directory is typically distributed as a .tar.gz archive that users
extract once before running TAtouScan.  Pass the extracted directory to
the ``--db`` CLI option.

The database version is checked at startup: databases older than
``MIN_DB_VERSION`` (or without ``version.txt``, i.e. released before 0.3.0)
are rejected.

Example::

    tatouscan --db /path/to/tatouscan_db/ --gff genome.gff --faa genome.faa …
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

# Oldest database version compatible with this version of TAtouScan.
MIN_DB_VERSION = "0.3.0"

# Zenodo concept record: always resolves to the latest database version.
DB_DOWNLOAD_URL = "https://doi.org/10.5281/zenodo.15305312"

VERSION_FILENAME = "version.txt"

# Fixed filenames inside the database directory.
DB_FILENAMES: dict[str, str] = {
    "hmm_db": "ta.hmm",
    "hmm_info": "hmm_info.tsv",
    "ref_stats": "family_statistics.tsv",
    "known_pairs": "known_pairs.tsv",
}


class DatabaseError(Exception):
    """Base class of the errors raised when a TAtouScan database cannot be used."""


class DatabaseVersionError(DatabaseError):
    """Raised when the database is missing a version or is too old."""


class DatabaseNotADirectoryError(DatabaseError, NotADirectoryError):
    """Raised when the database path does not exist or is not a directory."""


class DatabaseFileNotFoundError(DatabaseError, FileNotFoundError):
    """Raised when a required file is missing from the database directory."""


@dataclass(frozen=True)
class TAtouScanDB:
    """Resolved paths to the components of a TAtouScan database directory."""

    version: str
    hmm_db: Path
    hmm_info: Path
    ref_stats: Path
    known_pairs: Path


def parse_version(version: str) -> tuple[int, ...]:
    """Parse a ``MAJOR.MINOR.PATCH`` version string into a comparable tuple.

    Raises
    ------
    ValueError
        If *version* is not made of dot-separated integers.
    """
    return tuple(int(part) for part in version.strip().split("."))


def check_db_version(path: Path, min_version: str = MIN_DB_VERSION) -> str:
    """Read ``version.txt`` in *path* and check it is at least *min_version*.

    Returns the database version string.

    Raises
    ------
    DatabaseVersionError
        If the version file is missing, unreadable, or older than *min_version*.
    """
    update_hint = (
        f"Please download the latest TAtouScan database from {DB_DOWNLOAD_URL} "
        "and extract it with: tar -xzf tatouscan_db.tar.gz"
    )
    version_file = path / VERSION_FILENAME
    if not version_file.exists():
        raise DatabaseVersionError(
            f"TAtouScan database at '{path}' has no {VERSION_FILENAME} file: it was "
            f"released before version {min_version} and is not compatible with this "
            f"version of TAtouScan (requires database >= {min_version}).\n{update_hint}"
        )

    version = version_file.read_text().strip()
    try:
        parsed = parse_version(version)
    except ValueError:
        raise DatabaseVersionError(
            f"Invalid database version '{version}' in {version_file} "
            "(expected MAJOR.MINOR.PATCH).\n" + update_hint
        ) from None

    if parsed < parse_version(min_version):
        raise DatabaseVersionError(
            f"TAtouScan database at '{path}' is version {version}, but this version of "
            f"TAtouScan requires database >= {min_version}.\n{update_hint}"
        )
    return version


def load_db(path: Path) -> TAtouScanDB:
    """Validate *path* as a TAtouScan database directory and return resolved paths.

    Raises
    ------
    DatabaseNotADirectoryError
        If *path* does not exist or is not a directory.
    DatabaseVersionError
        If the database has no version file or is older than ``MIN_DB_VERSION``.
    DatabaseFileNotFoundError
        If any of the required files is missing from the directory.

    All these errors are subclasses of :class:`DatabaseError`.
    """
    if not path.exists():
        raise DatabaseNotADirectoryError(f"Database path does not exist: {path}")
    if not path.is_dir():
        raise DatabaseNotADirectoryError(
            f"--db must point to a directory (not a file): {path}\n"
            "If you have a .tar.gz archive, extract it first: tar -xzf <archive>"
        )

    version = check_db_version(path)

    missing = [fname for fname in DB_FILENAMES.values() if not (path / fname).exists()]
    if missing:
        raise DatabaseFileNotFoundError(
            f"TAtouScan database at '{path}' is missing required file(s): "
            + ", ".join(missing)
            + f"\nExpected: {', '.join(DB_FILENAMES.values())}"
        )

    return TAtouScanDB(
        version=version,
        hmm_db=path / DB_FILENAMES["hmm_db"],
        hmm_info=path / DB_FILENAMES["hmm_info"],
        ref_stats=path / DB_FILENAMES["ref_stats"],
        known_pairs=path / DB_FILENAMES["known_pairs"],
    )
