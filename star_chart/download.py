import configparser
import logging
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import geopandas as gpd
import pandas as pd
import typer
from shapely.geometry import LineString, Point, Polygon
from skyfield.api import Loader, Star, wgs84
from skyfield.data import hipparcos, stellarium
from skyfield.projections import build_stereographic_projection
from typing_extensions import Annotated

from star_chart.constants import ascii_to_greek_letters

logger = logging.getLogger(__name__)


def download(
    date: Annotated[datetime, typer.Argument()],
    latitude: Annotated[float, typer.Argument(min=-90, max=90)],
    longitude: Annotated[float, typer.Argument(min=-180, max=180)],
    config_file_path: Annotated[
        Optional[Path],
        typer.Option(
            "--config",
            "-c",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
        ),
    ] = Path("./config.ini"),
    download_folder_path: Annotated[
        Optional[Path],
        typer.Option(
            "--download",
            "-d",
            exists=False,
            file_okay=False,
            dir_okay=True,
            resolve_path=True,
        ),
    ] = Path("./download"),
    output_file_path: Annotated[
        Optional[Path],
        typer.Option(
            "--output",
            "-o",
            exists=False,
            file_okay=False,
            dir_okay=True,
            resolve_path=True,
        ),
    ] = Path("./output.gpkg"),
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
    logger.info("Options:")
    logger.info("  config_file_path: %s", config_file_path)
    logger.info("  download_folder_path: %s", download_folder_path)
    logger.info("  output_file_path: %s", output_file_path)

    # Read configuration file
    config_file = configparser.ConfigParser()
    config_file.read(config_file_path)

    # Get URLs from configuration file
    ephemeris_data_url = config_file["ephemeris"]["data_url"]

    star_data_url = config_file["star"]["data_url"]
    star_name_url = config_file["star"]["name_url"]

    constellation_boundary_url = config_file["constellation"]["boundary_url"]
    constellation_data_url = config_file["constellation"]["data_url"]
    constellation_name_url = config_file["constellation"]["name_url"]

    # Create download directory
    logger.info(f"Creating download directory at {download_folder_path}")
    os.makedirs(download_folder_path, exist_ok=True)

    # Create default file loader
    load = Loader(download_folder_path)

    # Load ephemeris data
    logger.info("Downloading ephemeris data")
    ephemeris = load(ephemeris_data_url)

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

    # Load Hipparcos data
    logger.info("Downloading Hipparcos data")
    with load.open(star_data_url) as downloaded_file:
        stars_df = hipparcos.load_dataframe(downloaded_file)

    # Compute the x and y coordinates based on the projection
    star_positions = earth.at(timescale).observe(Star.from_dataframe(stars_df))  # type: ignore
    stars_df["x"], stars_df["y"] = projection(star_positions)

    # Set global epoch year
    epoch_year = stars_df["epoch_year"].max()

    # Drop not used columns
    stars_df = stars_df.drop(
        columns=[
            "parallax_mas",
            "ra_mas_per_year",
            "dec_mas_per_year",
            "ra_hours",
            "epoch_year",
        ]
    )

    logger.info("Downloading HYG data")
    # Load HYG data (Hipparcos-Yale-Gliese) to get proper names, constellations and bayer
    with load.open(star_name_url) as downloaded_file:
        star_names_df = pd.read_csv(
            downloaded_file,
            index_col="hip",
            usecols=["hip", "proper", "bayer", "con"],
        )
    star_names_df.columns = ("proper", "bayer", "abbreviation")
    star_names_df["abbreviation"] = star_names_df["abbreviation"].str.upper()
    star_names_df["bayer"] = star_names_df["bayer"].str.split("-").str[0]
    star_names_df["greek"] = star_names_df["bayer"].map(ascii_to_greek_letters)

    # Join columns
    stars_df = stars_df.merge(
        star_names_df, how="inner", left_index=True, right_index=True
    )
    stars_df.index = stars_df.index.astype("int")

    # Build geodataframe
    logger.info("Build stars geodataframe")
    stars_gdf = gpd.GeoDataFrame(
        stars_df,
        geometry=gpd.points_from_xy(stars_df["x"], stars_df["y"]),
        crs=4326,
    )  # type: ignore
    stars_gdf = stars_gdf[
        ["magnitude", "proper", "bayer", "greek", "abbreviation", "geometry"]
    ]

    # Set CRS
    stars_gdf = stars_gdf.set_crs(epsg=4326)  # type: ignore

    # Export to .gpkg
    logger.info(f"Saving stars geodataframe to {output_file_path}")
    stars_gdf.to_file(output_file_path, driver="GPKG", layer="stars", encoding="utf-8")

    # Load constellation name data
    with load.open(constellation_name_url) as downloaded_file:
        constellations_names_df = pd.read_csv(
            downloaded_file,
            delimiter="\t",
            header=None,
            names=["abbreviation", "name"],
            index_col=0,
            usecols=[0, 1],
        )
    constellations_names_df.index = constellations_names_df.index.str.upper()

    # Load constellation edge data
    with load.open(constellation_data_url) as downloaded_file:
        constellations = stellarium.parse_constellations(downloaded_file)
    constellations_coordinate_pairs = sum(
        (
            [
                [constellation[0].upper(), vertice[0], vertice[1]]
                for vertice in constellation[1]
            ]
            for constellation in constellations
        ),
        [],
    )
    constellations_df = pd.DataFrame(
        constellations_coordinate_pairs,
        columns=["abbreviation", "start_hip", "end_hip"],
    )
    constellations_df.index.name = "id"

    # Merge with star coordinates
    constellations_df = constellations_df.merge(
        stars_df[["x", "y", "ra_degrees", "dec_degrees"]].add_prefix("start_"),
        how="inner",
        left_on="start_hip",
        right_index=True,
    )
    constellations_df = constellations_df.merge(
        stars_df[["x", "y", "ra_degrees", "dec_degrees"]].add_prefix("end_"),
        how="inner",
        left_on="end_hip",
        right_index=True,
    )

    # Build geodataframe
    constellations_gdf = gpd.GeoDataFrame(
        constellations_df,
        geometry=constellations_df.apply(
            lambda item: LineString(
                [
                    Point([item["start_x"], item["start_y"]]),
                    Point([item["end_x"], item["end_y"]]),
                ]
            ),  # type: ignore
            result_type=None,
            axis=1,
        ),  # type: ignore
        crs=4326,
    )  # type: ignore
    constellations_gdf = constellations_gdf.dissolve(by="abbreviation").merge(
        constellations_names_df, left_on="abbreviation", right_index=True
    )[["geometry", "name"]]

    # Set CRS
    constellations_gdf = constellations_gdf.set_crs(epsg=4326)  # type: ignore

    # Export to .gpkg
    constellations_gdf.to_file(
        output_file_path, driver="GPKG", layer="constellations", encoding="utf-8"
    )

    # Load constellation boundary data
    with load.open(constellation_boundary_url) as downloaded_file:
        colspecs = [(0, 9), (9, 18), (19, 23), (24, 25)]
        names = ["ra_hours", "dec_degrees", "abbreviation", "type_point"]
        boundaries_df = pd.read_fwf(
            downloaded_file,
            colspecs=colspecs,
            names=names,
        )

    # Filter boundary data by point type
    boundaries_df = boundaries_df[boundaries_df["type_point"] == "O"]

    # Create required columns
    boundaries_df.index.name = "id"
    boundaries_df["epoch_year"] = epoch_year

    # Convert hours to degrees
    boundaries_df["ra_degrees"] = boundaries_df["ra_hours"]

    # Project coordinates
    boundary_positions = earth.at(timescale).observe(Star.from_dataframe(boundaries_df))  # type: ignore
    boundaries_df["x"], boundaries_df["y"] = projection(boundary_positions)

    # Drop not used columns
    boundaries_df = boundaries_df.drop(columns=["type_point", "ra_hours", "epoch_year"])

    # Build geodataframe
    boundaries_gdf = gpd.GeoDataFrame(
        boundaries_df,
        geometry=gpd.points_from_xy(boundaries_df["x"], boundaries_df["y"]),
        crs=4326,
    )  # type: ignore
    boundaries_gdf = (
        boundaries_gdf.groupby("abbreviation")["geometry"]
        .apply(lambda x: Polygon(x.tolist()))
        .reset_index()
    )
    boundaries_gdf.index.name = "id"

    # Remove numbers from new abbreviations
    boundaries_gdf["new_abbreviation"] = boundaries_gdf["abbreviation"].str[:3]

    # Add constellation names
    boundaries_gdf = boundaries_gdf.merge(
        constellations_names_df, left_on="new_abbreviation", right_index=True
    )[["geometry", "name", "new_abbreviation", "abbreviation"]]
    boundaries_gdf = boundaries_gdf.drop(columns=["new_abbreviation"])
    boundaries_gdf = boundaries_gdf.set_index("abbreviation", verify_integrity=True)

    # Set CRS
    boundaries_gdf = boundaries_gdf.set_crs(epsg=4326)  # type: ignore

    # Export to .gpkg
    boundaries_gdf.to_file(
        output_file_path, driver="GPKG", layer="boundaries", encoding="utf-8"
    )
