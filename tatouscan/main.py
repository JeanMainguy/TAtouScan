from typing import Optional

import typer
from typing_extensions import Annotated

from tatouscan import __version__
from rich.logging import RichHandler
import logging

from tatouscan.arguments import validate_args, args, logging_levels

logger = logging.getLogger(__name__)

app = typer.Typer(
    help="TAtouScan: A tool for identifying toxin-antitoxin (TA) systems."
)


def setup_logging(verbose_mode: str):
    """Setup logging for TAtouScan"""
    logging_level = logging_levels[verbose_mode]
    logging.basicConfig(
        level=logging_level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler()],
    )


def show_wip_message():
    """Print 'under development' warning"""
    typer.echo(
        "TAtouScan CLI is under development. Run `tatouscan --help` for available commands.",
        color=True,
    )


@app.command(no_args_is_help=True)
def main(
    faa_file: Annotated[str, args["faa_file"]],
    gff_file: Annotated[str, args["gff_file"]],
    hmm_db: Annotated[str, args["hmm_db"]],
    output_directory: Annotated[str, args["output_directory"]],
    e_value: Annotated[float, args["e_value"]],
    max_sequence_length: Annotated[int, args["max_sequence_length"]],
    max_distance: Annotated[int, args["max_distance"]],
    cpu: Annotated[int, args["cpu"]] = 1,
    verbose_mode: Annotated[str, args["verbose_mode"]] = "INFO",
    version: Annotated[Optional[bool], args["version"]] = None,
):
    """Main entry point for TAtouScan CLI."""
    validate_args(faa_file, gff_file, hmm_db, output_directory, e_value, max_sequence_length, max_distance, cpu, verbose_mode)
    setup_logging(verbose_mode)
    show_wip_message()
    tatouscan(faa_file, gff_file, hmm_db, output_directory, e_value, max_sequence_length, max_distance, cpu)


if __name__ == "__main__":
    app()
