from pathlib import Path

import pytest

from tatouscan.database import (
    DB_FILENAMES,
    MIN_DB_VERSION,
    VERSION_FILENAME,
    DatabaseVersionError,
    check_db_version,
    load_db,
    parse_version,
)


def make_db(tmp_path: Path, version: str | None = MIN_DB_VERSION) -> Path:
    """Create a minimal database directory with empty files."""
    db = tmp_path / "tatouscan_db"
    db.mkdir()
    for fname in DB_FILENAMES.values():
        (db / fname).write_text("")
    if version is not None:
        (db / VERSION_FILENAME).write_text(f"{version}\n")
    return db


def test_parse_version():
    assert parse_version("0.3.0") == (0, 3, 0)
    assert parse_version(" 1.10.2\n") == (1, 10, 2)
    assert parse_version("0.10.0") > parse_version("0.9.0")
    with pytest.raises(ValueError):
        parse_version("v0.3")


def test_load_db_current_version(tmp_path: Path):
    db = load_db(make_db(tmp_path))
    assert db.version == MIN_DB_VERSION
    assert db.hmm_db.name == "ta.hmm"


def test_load_db_newer_version(tmp_path: Path):
    major, minor, patch = parse_version(MIN_DB_VERSION)
    newer = f"{major}.{minor}.{patch + 1}"
    assert load_db(make_db(tmp_path, newer)).version == newer


def test_load_db_without_version_file(tmp_path: Path):
    """Databases released before 0.3.0 have no version.txt and are rejected."""
    with pytest.raises(DatabaseVersionError, match="no version.txt"):
        load_db(make_db(tmp_path, version=None))


def test_check_db_version_too_old(tmp_path: Path):
    db = make_db(tmp_path, "0.2.0")
    with pytest.raises(DatabaseVersionError, match="requires database >= 0.3.0"):
        check_db_version(db, min_version="0.3.0")


def test_check_db_version_invalid(tmp_path: Path):
    with pytest.raises(DatabaseVersionError, match="Invalid database version"):
        check_db_version(make_db(tmp_path, "latest"))


def test_load_db_missing_file(tmp_path: Path):
    db = make_db(tmp_path)
    (db / "known_pairs.tsv").unlink()
    with pytest.raises(FileNotFoundError, match="known_pairs.tsv"):
        load_db(db)
