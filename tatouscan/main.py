from typing import Optional

import typer
from typing_extensions import Annotated

from tatouscan import __version__
from rich.logging import RichHandler
import logging
from pathlib import Path

from tatouscan.annotation import annotate_cdss

logger = logging.getLogger(__name__)

app = typer.Typer(
    help="TAtouScan: A tool for identifying toxin-antitoxin (TA) systems."
)


def version_callback(value: bool):
    """Prints the version and exits if --version is passed."""
    if value:
        typer.echo(f"TAtouScan {__version__}")
        raise typer.Exit()


@app.command(no_args_is_help=True)
def main(
    gff: Path,
    faa: Path,
    hmm_db: Path,
    version: Annotated[
        Optional[bool],
        typer.Option(
            "--version", callback=version_callback, help="Show the version and exit."
        ),
    ] = None,
):

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler()],
    )

    """Main entry point for TAtouScan CLI."""
    typer.echo(
        "TAtouScan CLI is under development. Run `tatouscan --help` for available commands.",
        color=True,
    )

    contig_name_and_cdss = annotate_cdss(
        gff_file=gff, faa_file=faa, hmm_db=hmm_db, e_value_threshold=0.01
    )

    for contig_name, cdss in contig_name_and_cdss:
        print(contig_name, len(cdss))
        for cds in cdss:
            if len(cds.ta_hits) > 0:
                print(cds.id, cds.coordinates, len(cds.ta_hits))

        print("=" * 80)


if __name__ == "__main__":
    app()
