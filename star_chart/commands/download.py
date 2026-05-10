import logging
from pathlib import Path
from typing import Any, Generator

import duckdb
import pandas as pd
import requests
import typer
from astroquery.vizier import Vizier
from astroquery.vizier.core import VizierClass
from pandas import DataFrame
from skyfield.api import load
from typing_extensions import Annotated

from star_chart.constants import (
    CONSTELLATION_EDGES_DOWNLOAD_URL,
    EPHEMERIS_DOWNLOAD_URL,
    EPOCH_YEAR,
)

logger = logging.getLogger(__name__)


def _get_stars_data_frame(vizier: VizierClass) -> DataFrame:
    DATASET = "I/239/hip_main"

    results = vizier.query_constraints(  # type: ignore
        catalog=DATASET,
        **{"_RA.icrs": "!= NULL", "_DE.icrs": "!= NULL"},
    )
    result = next(iter(results), None)
    if result is None:
        raise RuntimeError(f"{DATASET} not found.")

    star_column_lookup = {
        "HIP": "hip",
        "Vmag": "magnitude",
        "Plx": "parallax_mas",
        "pmRA": "ra_mas_per_year",
        "pmDE": "dec_mas_per_year",
        "_RA.icrs": "ra_degrees",
        "_DE.icrs": "dec_degrees",
    }

    filtered_result = result[list(star_column_lookup.keys())]
    for old_name, new_name in star_column_lookup.items():
        filtered_result.rename_column(old_name, new_name)
    filtered_result["ra_hours"] = filtered_result["ra_degrees"] / 15.0
    filtered_result.remove_column("ra_degrees")
    filtered_result["epoch_year"] = EPOCH_YEAR

    stars_df = filtered_result.to_pandas(index="hip")
    return stars_df.reset_index()


def _get_constellation_boundaries_data_frame(vizier: VizierClass) -> DataFrame:
    DATASET = "VI/49/bound_18"

    constellation_boundary_column_lookup = {
        "cst": "iau_abbreviation",
        "_RA.icrs": "ra_degrees",
        "_DE.icrs": "dec_degrees",
    }

    results = vizier.query_constraints(  # type: ignore
        catalog=DATASET,
        **{"type": "O", "_RA.icrs": "!= NULL", "_DE.icrs": "!= NULL"},
    )
    result = next(iter(results), None)
    if result is None:
        raise RuntimeError(
            "Invalid Vizier output columns and/or constraints on columns."
        )

    filtered_result = result[list(constellation_boundary_column_lookup.keys())]
    for old_name, new_name in constellation_boundary_column_lookup.items():
        filtered_result.rename_column(old_name, new_name)
    filtered_result["ra_hours"] = filtered_result["ra_degrees"] / 15.0
    filtered_result.remove_columns(["ra_degrees"])
    filtered_result_df = filtered_result.to_pandas()

    # Duplicating first vertex of each constellation to close the boundary polygon
    constellation_boundaries_df = (
        filtered_result_df.groupby("iau_abbreviation", sort=False)
        .apply(lambda group: pd.concat([group, group.head(1)]), include_groups=False)
        .droplevel(1)
        .reset_index()
    )

    constellation_boundaries_df = constellation_boundaries_df.assign(
        order=constellation_boundaries_df.groupby("iau_abbreviation").cumcount()
    )

    constellation_boundaries_df["epoch_year"] = EPOCH_YEAR
    constellation_boundaries_df.index = constellation_boundaries_df.index.rename("id")
    constellation_boundaries_df = constellation_boundaries_df.reset_index()
    constellation_boundaries_df = constellation_boundaries_df.astype(
        {"id": "int32", "order": "int32"}
    )
    return constellation_boundaries_df


def _get_constellation_edges_data_frame():
    def _get_constellation_edges() -> Generator[tuple[str, str, int, int], Any, None]:
        response = requests.get(CONSTELLATION_EDGES_DOWNLOAD_URL).json()
        constellations = response["constellations"]
        for constellation in constellations:
            name = str(constellation["common_name"]["english"])
            abbreviation = str(constellation["iau"]).upper()
            for segment, line in enumerate(constellation["lines"]):
                for star in line:
                    if not isinstance(star, int):
                        continue
                    yield (abbreviation, name, segment, star)

    constellation_edges_df = pd.DataFrame.from_records(
        _get_constellation_edges(),
        columns=["iau_abbreviation", "iau_name", "segment", "hip"],
    )

    constellation_edges_df["epoch_year"] = EPOCH_YEAR
    constellation_edges_df.index = constellation_edges_df.index.rename("id")
    constellation_edges_df = constellation_edges_df.reset_index()
    constellation_edges_df = constellation_edges_df.astype(
        {"id": "int32", "segment": "int32", "hip": "int32"}
    )
    return constellation_edges_df


def _download_ephemeris(output_ephemeris_file_path: Path):
    if not output_ephemeris_file_path.exists():
        load.download(EPHEMERIS_DOWNLOAD_URL, str(output_ephemeris_file_path))


def download(
    output_database_file_path: Annotated[
        Path,
        typer.Argument(
            help="Output DuckDB database file path.",
            file_okay=True,
            dir_okay=False,
            resolve_path=True,
        ),
    ] = Path("./output/download.duckdb"),
    output_ephemeris_file_path: Annotated[
        Path,
        typer.Argument(
            help="Output ephemeris file path.",
            file_okay=True,
            dir_okay=False,
            resolve_path=True,
        ),
    ] = Path("./output/de421.bsp"),
):
    """Automatically download, sanitize, and organize astronomical observation data from
    VizieR (https://vizier.cds.unistra.fr/) and Stellarium (https://github.com/Stellarium)
    into a DuckDB database file containing the following tables: stars,
    constellation_boundaries and constellation_edges. Ephemeris data from NASA's JPL
    will also be downloaded into a BSP file.
    """

    logger.info("Arguments:")
    logger.info("  output_database_file_path: %s", output_database_file_path)
    logger.info("  output_ephemeris_file_path: %s", output_ephemeris_file_path)

    vizier = Vizier()
    vizier.ROW_LIMIT = -1

    logger.info("Downloading stars")
    stars_df = _get_stars_data_frame(vizier)  # noqa: F841
    logger.info("Downloading constellation boundaries")
    constellation_boundaries_df = _get_constellation_boundaries_data_frame(vizier)  # noqa: F841
    logger.info("Downloading constellation edges")
    constellation_edges_df = _get_constellation_edges_data_frame()  # noqa: F841

    with duckdb.connect(output_database_file_path) as connection:
        logger.info("Creating 'stars' table")
        connection.execute(
            "CREATE OR REPLACE TABLE stars AS SELECT * FROM stars_df",
        )
        logger.info("Creating 'constellation_boundaries' table")
        connection.execute(
            "CREATE OR REPLACE TABLE constellation_boundaries AS SELECT * FROM constellation_boundaries_df",
        )
        logger.info("Creating 'constellation_edges' table")
        connection.execute(
            "CREATE OR REPLACE TABLE constellation_edges AS SELECT * FROM constellation_edges_df",
        )

    _download_ephemeris(output_ephemeris_file_path)
