import logging
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import geopandas as gpd
import typer
from geopandas import GeoDataFrame
from pandas import DataFrame
from skyfield.api import Star, load, wgs84
from skyfield.projections import build_stereographic_projection
from typing_extensions import Annotated

logger = logging.getLogger(__name__)

# Set global epoch year (1991.25)
EPOCH_YEAR = 1991.25


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
        stars_df,
        geometry=gpd.points_from_xy(stars_df["x"], stars_df["y"]),
        crs=4326,
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
    return


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
    return


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

    stars_gdf = _build_stars_geo_data_frame(
        earth,
        timescale,
        projection,
        stars_df,
    )
    constellation_edges_gdf = _build_constellation_edges_geo_data_frame(
        earth,
        timescale,
        projection,
        stars_df,
        constellation_edges_df,
    )

    constellation_boundaries_gdf = _build_constellation_boundaries_geo_data_frame(
        earth,
        timescale,
        projection,
        constellation_edges_df,
        constellation_boundaries_df,
    )

    # logger.info("Downloading HYG data")
    # # Load HYG data (Hipparcos-Yale-Gliese) to get proper names, constellations and bayer
    # with load.open(star_name_url) as compressed_downloaded_file:
    #     extracted_downloaded_file = extract_file(compressed_downloaded_file)
    #     star_names_df = pd.read_csv(
    #         extracted_downloaded_file,
    #         index_col="hip",
    #         usecols=["hip", "proper", "bayer", "con"],
    #     )
    # star_names_df.columns = ("proper", "bayer", "abbreviation")
    # star_names_df["abbreviation"] = star_names_df["abbreviation"].str.upper()
    # star_names_df["bayer"] = star_names_df["bayer"].str.split("-").str[0]
    # star_names_df["greek"] = star_names_df["bayer"].map(ascii_to_greek_letters)

    # # Join columns
    # stars_df = stars_df.merge(
    #     star_names_df, how="inner", left_index=True, right_index=True
    # )
    # stars_df.index = stars_df.index.astype("int")

    # # Build geodataframe
    # logger.info("Building stars geodataframe")
    # stars_gdf = gpd.GeoDataFrame(
    #     stars_df,
    #     geometry=gpd.points_from_xy(stars_df["x"], stars_df["y"]),
    #     crs=4326,
    # )  # type: ignore
    # stars_gdf = stars_gdf[
    #     ["magnitude", "proper", "bayer", "greek", "abbreviation", "geometry"]
    # ]

    # # Set CRS
    # stars_gdf = stars_gdf.set_crs(epsg=4326)  # type: ignore

    # # Export to .gpkg
    # logger.info(f"Saving stars geodataframe to {output_file_path}")
    # stars_gdf.to_file(output_file_path, driver="GPKG", layer="stars", encoding="utf-8")

    # # Load constellation name data
    # logger.info("Downloading constellation names")
    # with load.open(constellation_name_url) as downloaded_file:
    #     constellations_names_df = pd.read_csv(
    #         downloaded_file,
    #         delimiter="\t",
    #         header=None,
    #         names=["abbreviation", "name"],
    #         index_col=0,
    #         usecols=[0, 1],
    #     )
    # constellations_names_df.index = constellations_names_df.index.str.upper()

    # # Load constellation edge data
    # logger.info("Downloading constellation edges")
    # with load.open(constellation_data_url) as downloaded_file:
    #     constellations = stellarium.parse_constellations(downloaded_file)
    # constellations_coordinate_pairs = sum(
    #     (
    #         [
    #             [constellation[0].upper(), vertice[0], vertice[1]]
    #             for vertice in constellation[1]
    #         ]
    #         for constellation in constellations
    #     ),
    #     [],
    # )
    # constellations_df = pd.DataFrame(
    #     constellations_coordinate_pairs,
    #     columns=["abbreviation", "start_hip", "end_hip"],
    # )
    # constellations_df.index.name = "id"

    # # Merge with star coordinates
    # logger.info("Obtaining constellation edge coordinates from star's coordinates")
    # constellations_df = constellations_df.merge(
    #     stars_df[["x", "y", "ra_degrees", "dec_degrees"]].add_prefix("start_"),
    #     how="inner",
    #     left_on="start_hip",
    #     right_index=True,
    # )
    # constellations_df = constellations_df.merge(
    #     stars_df[["x", "y", "ra_degrees", "dec_degrees"]].add_prefix("end_"),
    #     how="inner",
    #     left_on="end_hip",
    #     right_index=True,
    # )

    # # Build geodataframe
    # logger.info("Building constellation edges geodataframe")
    # constellations_gdf = gpd.GeoDataFrame(
    #     constellations_df,
    #     geometry=constellations_df.apply(
    #         lambda item: LineString(
    #             [
    #                 Point([item["start_x"], item["start_y"]]),
    #                 Point([item["end_x"], item["end_y"]]),
    #             ]
    #         ),  # type: ignore
    #         result_type=None,
    #         axis=1,
    #     ),  # type: ignore
    #     crs=4326,
    # )  # type: ignore
    # constellations_gdf = constellations_gdf.dissolve(by="abbreviation").merge(
    #     constellations_names_df, left_on="abbreviation", right_index=True
    # )[["geometry", "name"]]

    # # Set CRS
    # constellations_gdf = constellations_gdf.set_crs(epsg=4326)  # type: ignore

    # # Export to .gpkg
    # logger.info(f"Saving constellation eges geodataframe to {output_file_path}")
    # constellations_gdf.to_file(
    #     output_file_path, driver="GPKG", layer="constellations", encoding="utf-8"
    # )

    # # Load constellation boundary data
    # logger.info("Downloading constellation boundaries")
    # with load.open(constellation_boundary_url) as downloaded_file:
    #     colspecs = [(0, 9), (9, 18), (19, 23), (24, 25)]
    #     names = ["ra_hours", "dec_degrees", "abbreviation", "type_point"]
    #     boundaries_df = pd.read_fwf(
    #         downloaded_file,
    #         colspecs=colspecs,
    #         names=names,
    #     )

    # # Filter boundary data by point type
    # boundaries_df = boundaries_df[boundaries_df["type_point"] == "O"]

    # # Create required columns
    # boundaries_df.index.name = "id"
    # boundaries_df["epoch_year"] = EPOCH_YEAR

    # # Convert hours to degrees
    # boundaries_df["ra_degrees"] = boundaries_df["ra_hours"]

    # # Project coordinates
    # logger.info("Projecting constellation boundaries coordinates")
    # boundary_positions = earth.at(timescale).observe(Star.from_dataframe(boundaries_df))  # type: ignore
    # boundaries_df["x"], boundaries_df["y"] = projection(boundary_positions)

    # # Drop not used columns
    # boundaries_df = boundaries_df.drop(columns=["type_point", "ra_hours", "epoch_year"])

    # # Build geodataframe
    # logger.info("Building constellation boundaries geodataframe")
    # boundaries_gdf = gpd.GeoDataFrame(
    #     boundaries_df,
    #     geometry=gpd.points_from_xy(boundaries_df["x"], boundaries_df["y"]),
    #     crs=4326,
    # )  # type: ignore
    # boundaries_gdf = (
    #     boundaries_gdf.groupby("abbreviation")["geometry"]
    #     .apply(lambda x: Polygon(x.tolist()))
    #     .reset_index()
    # )
    # boundaries_gdf.index.name = "id"

    # # Remove numbers from new abbreviations
    # boundaries_gdf["new_abbreviation"] = boundaries_gdf["abbreviation"].str[:3]

    # # Add constellation names
    # boundaries_gdf = boundaries_gdf.merge(
    #     constellations_names_df, left_on="new_abbreviation", right_index=True
    # )[["geometry", "name", "new_abbreviation", "abbreviation"]]
    # boundaries_gdf = boundaries_gdf.drop(columns=["new_abbreviation"])
    # boundaries_gdf = boundaries_gdf.set_index("abbreviation", verify_integrity=True)

    # # Set CRS
    # boundaries_gdf = boundaries_gdf.set_crs(epsg=4326)  # type: ignore

    # # Export to .gpkg
    # logger.info(f"Saving constellation boundaries geodataframe to {output_file_path}")
    # boundaries_gdf.to_file(
    #     output_file_path, driver="GPKG", layer="boundaries", encoding="utf-8"
    # )
