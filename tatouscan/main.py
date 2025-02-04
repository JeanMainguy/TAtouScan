from typing import Optional

import typer
from typing_extensions import Annotated

from tatouscan import __version__

app = typer.Typer(
    help="TatouScan: A tool for identifying toxin-antitoxin (TA) systems."
)


def version_callback(value: bool):
    """Prints the version and exits if --version is passed."""
    if value:
        typer.echo(f"TatouScan {__version__}")
        raise typer.Exit()


@app.command()
def main(
    version: Annotated[
        Optional[bool],
        typer.Option(
            "--version", callback=version_callback, help="Show the version and exit."
        ),
    ] = None,
):
    """Main entry point for TatouScan CLI."""
    typer.echo(
        "TatouScan CLI is under development. Run `tatouscan --help` for available commands."
    )


if __name__ == "__main__":
    app()
