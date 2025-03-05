import logging
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import geopandas as gpd
import typer
from geopandas import GeoDataFrame
from pandas import DataFrame
from shapely import LineString, MultiLineString, Polygon
from skyfield.api import Star, load, wgs84
from skyfield.projections import build_stereographic_projection
from typing_extensions import Annotated

logger = logging.getLogger(__name__)

# Set global epoch year (1991.25)
EPOCH_YEAR = 1991.25

# Default EPSG code
EPSG_CODE = 3857


def _build_stars_geo_data_frame(
    earth,
    timescale,
    projection,
    stars_df: DataFrame,
) -> GeoDataFrame:
    # Compute the stars x and y coordinates based on the projection
    logger.info("Projecting 'stars' coordinates")
    star_positions = earth.at(timescale).observe(Star.from_dataframe(stars_df))  # type: ignore
    stars_df["x"], stars_df["y"] = projection(star_positions)

    logger.info("Building stars geodataframe")
    stars_gdf = gpd.GeoDataFrame(
        stars_df[["hip", "magnitude"]].set_index("hip"),
        geometry=gpd.points_from_xy(stars_df["x"], stars_df["y"]),
        crs=EPSG_CODE,
    )
    return stars_gdf


def _build_constellation_edges_geo_data_frame(
    earth,
    timescale,
    projection,
    stars_df: DataFrame,
    constellation_edges_df: DataFrame,
) -> GeoDataFrame:
    logger.info("Merging 'constellation_edges' with the 'stars' table for coordinates")
    constellation_edges_df = constellation_edges_df.merge(
        stars_df[
            [
                "hip",
                "dec_degrees",
                "ra_hours",
                "parallax_mas",
                "ra_mas_per_year",
                "dec_mas_per_year",
            ]
        ],
        how="inner",
        left_on="hip",
        right_on="hip",
    )

    # Compute the constellation edges x and y coordinates based on the projection
    logger.info("Projecting 'constellation_edges' coordinates")
    constellation_edges_positions = earth.at(timescale).observe(
        Star.from_dataframe(constellation_edges_df)
    )  # type: ignore
    constellation_edges_df["x"], constellation_edges_df["y"] = projection(
        constellation_edges_positions
    )

    constellation_edges_grouped_series = (
        constellation_edges_df.groupby(["iau_abbreviation", "iau_name", "segment"])
        .apply(lambda group: LineString(zip(group["x"], group["y"])))
        .groupby(["iau_abbreviation", "iau_name"])
        .apply(lambda group: MultiLineString(list(group)))
    )

    constellation_edges_grouped_df = constellation_edges_grouped_series.to_frame(
        name="geometry"
    )
    constellation_edges_gdf = GeoDataFrame(
        constellation_edges_grouped_df,
        geometry="geometry",
        crs=EPSG_CODE,
    )
    return constellation_edges_gdf


def _build_constellation_boundaries_geo_data_frame(
    earth,
    timescale,
    projection,
    constellation_edges_df: DataFrame,
    constellation_boundaries_df: DataFrame,
) -> GeoDataFrame:
    logger.info("Merging 'constellation_edges' with the 'stars' table for coordinates")
    constellation_boundaries_df = constellation_boundaries_df.merge(
        constellation_edges_df[
            [
                "iau_abbreviation",
                "iau_name",
            ]
        ],
        how="inner",
        left_on="iau_abbreviation",
        right_on="iau_abbreviation",
    )

    # Compute the constellation boundaries x and y coordinates based on the projection
    logger.info("Projecting 'constellation_boundaries' coordinates")
    constellation_boundaries_positions = earth.at(timescale).observe(
        Star.from_dataframe(constellation_boundaries_df)
    )  # type: ignore
    constellation_boundaries_df["x"], constellation_boundaries_df["y"] = projection(
        constellation_boundaries_positions
    )

    constellation_boundaries_grouped_series = constellation_boundaries_df.groupby(
        ["iau_abbreviation", "iau_name"]
    ).apply(lambda group: Polygon(list(zip(group["x"], group["y"]))))

    constellation_boundaries_grouped_df = (
        constellation_boundaries_grouped_series.to_frame(name="geometry")
    )
    constellation_boundaries_gdf = GeoDataFrame(
        constellation_boundaries_grouped_df,
        geometry="geometry",
        crs=EPSG_CODE,
    )

    return constellation_boundaries_gdf


def observe(
    date: Annotated[datetime, typer.Argument()],
    latitude: Annotated[float, typer.Argument(min=-90, max=90)],
    longitude: Annotated[float, typer.Argument(min=-180, max=180)],
    input_database_file_path: Annotated[
        Path,
        typer.Argument(
            help="Input DuckDB database file path.",
            file_okay=True,
            dir_okay=False,
            resolve_path=True,
        ),
    ] = Path("./download.duckdb"),
    input_ephemeris_file_path: Annotated[
        Path,
        typer.Argument(
            help="Input ephemeris file path.",
            file_okay=True,
            dir_okay=False,
            resolve_path=True,
        ),
    ] = Path("./de421.bsp"),
    output_geopackage_file_path: Annotated[
        Path,
        typer.Argument(
            help="Output GeoPackage file path.",
            file_okay=True,
            dir_okay=False,
            resolve_path=True,
        ),
    ] = Path("./download.gpkg"),
):
    """Automatically downloads, sanitizes, and organizes astronomical observation data
    from a given DATE and location (LATITUDE, LONGITUDE) into a GeoPackage. Data
    providers are defined in the CONFIG file.

    For convenience, the specified temporary download directory option will be created
    if nonexistent.
    """

    logger.info("Arguments:")
    logger.info("  date: %s", date)
    logger.info("  latitude: %s", latitude)
    logger.info("  longitude: %s", longitude)
    logger.info("  input database: %s", input_database_file_path)
    logger.info("  ephemeris: %s", input_ephemeris_file_path)
    logger.info("  output geopackage: %s", output_geopackage_file_path)

    # Load ephemeris data
    logger.info(f"Loading ephemeris from {input_database_file_path}")
    ephemeris = load(str(input_ephemeris_file_path))

    # An ephemeris on Earth position
    earth = ephemeris["earth"]

    # Get specified date as an UTC-based datetime object
    utc = date.replace(tzinfo=timezone.utc)

    # Define observer using specified location coordinates and UTC time
    timescale = load.timescale().from_datetime(utc)
    observer = wgs84.latlon(latitude, longitude).at(timescale)

    # Define a center based on observer position
    ra, dec, _ = observer.radec()
    center = earth.at(timescale).observe(Star(ra=ra, dec=dec))  # type: ignore

    # Build the stereographic projection
    logger.info("Building the stereographic projection")
    projection = build_stereographic_projection(center)

    logger.info(
        "Connecting to the input DuckDB database at %s", input_database_file_path
    )
    with duckdb.connect(input_database_file_path) as input_database_connection:
        # Load star data
        logger.info("Loading 'stars' table")
        stars_df = input_database_connection.sql(
            "SELECT * FROM stars",
        ).df()

        # Load constellation edges data
        logger.info("Loading 'constellation_edges' table")
        constellation_edges_df = input_database_connection.sql(
            "SELECT * FROM constellation_edges"
        ).df()

        # Load constellation boundaries data
        logger.info("Loading 'constellation_boundaries' table")
        constellation_boundaries_df = input_database_connection.sql(
            "SELECT * FROM constellation_boundaries"
        ).df()

    logger.info("Building 'stars' GeoDataFrame")
    stars_gdf = _build_stars_geo_data_frame(
        earth,
        timescale,
        projection,
        stars_df,
    )
    logger.info("Saving 'stars' GeoDataFrame to %s/stars", output_geopackage_file_path)
    stars_gdf.to_file(
        output_geopackage_file_path,
        driver="GPKG",
        layer="stars",
        encoding="utf-8",
    )

    logger.info("Building 'constellation_edges' GeoDataFrame")
    constellation_edges_gdf = _build_constellation_edges_geo_data_frame(
        earth,
        timescale,
        projection,
        stars_df,
        constellation_edges_df,
    )
    logger.info(
        "Saving 'constellation_edges' GeoDataFrame to %s/constellation_edges",
        output_geopackage_file_path,
    )
    constellation_edges_gdf.to_file(
        output_geopackage_file_path,
        driver="GPKG",
        layer="constellation_edges",
        encoding="utf-8",
    )

    logger.info("Building 'constellation_boundaries' GeoDataFrame")
    constellation_boundaries_gdf = _build_constellation_boundaries_geo_data_frame(
        earth,
        timescale,
        projection,
        constellation_edges_df,
        constellation_boundaries_df,
    )
    logger.info(
        "Saving 'constellation_boundaries' GeoDataFrame to %s/constellation_boundaries",
        output_geopackage_file_path,
    )
    constellation_boundaries_gdf.to_file(
        output_geopackage_file_path,
        driver="GPKG",
        layer="constellation_boundaries",
        encoding="utf-8",
    )
