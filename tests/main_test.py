from pathlib import Path

import pytest
from typer.testing import CliRunner

from tatouscan.database import DB_FILENAMES, MIN_DB_VERSION, VERSION_FILENAME
from tatouscan.main import app

runner = CliRunner()


def test_cli_reports_missing_database_file_without_traceback(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
):
    db = tmp_path / "tatouscan_db"
    db.mkdir()
    (db / VERSION_FILENAME).write_text(f"{MIN_DB_VERSION}\n")
    for fname in DB_FILENAMES.values():
        if fname != "known_pairs.tsv":
            (db / fname).write_text("")

    gff, faa = tmp_path / "genome.gff", tmp_path / "genome.faa"
    gff.write_text("##gff-version 3\n")
    faa.write_text(">protein1\nMTEITAAMVK\n")

    result = runner.invoke(
        app,
        [
            "--gff",
            str(gff),
            "--faa",
            str(faa),
            "--db",
            str(db),
            "-o",
            str(tmp_path / "out"),
        ],
    )

    assert result.exit_code == 1
    # Handled error: the CLI exits cleanly instead of raising the exception
    assert result.exception is None or isinstance(result.exception, SystemExit)
    # The error is logged, not printed to stdout
    assert "missing required file(s): known_pairs.tsv" in caplog.text
