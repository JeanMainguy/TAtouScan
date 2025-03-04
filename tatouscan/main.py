from typing import Optional

import typer
from typing_extensions import Annotated

from tatouscan import __version__
from rich.logging import RichHandler
import logging
from pathlib import Path

from tatouscan.parser import parse_gff_file

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

    contig_to_cds = parse_gff_file(gff)
    contig, cds_list = list(contig_to_cds).pop()

    print(contig)
    for cds in cds_list:
        print(cds)


if __name__ == "__main__":
    app()
