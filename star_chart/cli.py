import logging

import typer

from star_chart.download import download
from star_chart.observe import observe
from star_chart.plot import plot

logging.basicConfig(
    format="%(asctime)s|%(levelname)8s| %(message)s",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

app = typer.Typer()
app.command()(download)
app.command()(observe)
app.command()(plot)

if __name__ == "__main__":
    app()
