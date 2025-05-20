from typing import Optional

import typer
from typing_extensions import Annotated

from tatouscan import __version__
from rich.logging import RichHandler
import logging
from pathlib import Path

from tatouscan.annotation import annotate_cdss, parse_hmm_db_info
from tatouscan.system import group_cdss_with_ta_annotation
from tatouscan.writer import write_gene_with_ta_annotation

logger = logging.getLogger(__name__)

app = typer.Typer(
    help="TAtouScan: A tool for identifying toxin-antitoxin (TA) systems.",
    context_settings={"help_option_names": ["-h", "--help"]},
)


def version_callback(value: bool):
    """Prints the version and exits if --version is passed."""
    if value:
        typer.echo(f"TAtouScan {__version__}")
        raise typer.Exit()


@app.command(no_args_is_help=True)
def main(
    gff: Annotated[
        Path,
        typer.Option(
            "--gff",
            help="Path to the GFF file containing gene annotations.",
            exists=True,
        ),
    ],
    faa: Annotated[
        Path,
        typer.Option(
            "--faa",
            help="Path to the FASTA file containing protein sequences.",
            exists=True,
        ),
    ],
    hmm_db: Annotated[
        Path,
        typer.Option(
            "--hmm_db",
            help="Path to the HMM database file.",
            exists=True,
        ),
    ],
    hmm_info: Annotated[
        Path,
        typer.Option(
            "--hmm_info",
            help="Path to a TSV containing informaiton on the HMM profile.",
            exists=True,
        ),
    ],
    output_file: Annotated[
        Path,
        typer.Option(
            "--output",
            help="Path to a TSV file to write the gene with TA annotation.",
        ),
    ] = Path("tatouscan_results.tsv"),
    max_distance: Annotated[
        int,
        typer.Option(
            help="Maximum distance in nucleotides between genes to be considered in the same TA cluster."
        ),
    ] = 500,
    max_e_value: Annotated[
        float,
        typer.Option(
            help="Maximum E-value for TA hits to be considered.",
        ),
    ] = 0.005,
    version: Annotated[
        Optional[bool],
        typer.Option(
            "--version",
            callback=version_callback,
            help="Show the version and exit.",
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
    cds_protein_attr = ["id", "protein_id", "name", "locus_tag"]

    contig_name_and_cdss = annotate_cdss(
        gff_file=gff,
        faa_file=faa,
        hmm_db=hmm_db,
        e_value_threshold=max_e_value,
        matching_attributes=cds_protein_attr,
    )

    contig_name_and_cdss_with_ta_hit = group_cdss_with_ta_annotation(
        contig_name_and_cdss, max_distance
    )

    hmm_name_to_info = parse_hmm_db_info(hmm_info_file=hmm_info)

    write_gene_with_ta_annotation(
        contig_name_and_cdss_with_ta_hit, hmm_name_to_info, output_file
    )


if __name__ == "__main__":
    app()
