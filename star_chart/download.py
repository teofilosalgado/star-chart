import pandas as pd
from astroquery.vizier import Vizier
from astroquery.vizier.core import VizierClass
from pandas import DataFrame

EPOCH_YEAR = 1991.25


def _get_stars_data_frame(vizier: VizierClass) -> DataFrame:
    DATASET = "I/239/hip_main"

    results = vizier.query_constraints(
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
        "_DE.icrs": "de_degrees",
    }

    filtered_result = result[list(star_column_lookup.keys())]
    for old_name, new_name in star_column_lookup.items():
        filtered_result.rename_column(old_name, new_name)
    filtered_result["ra_hours"] = filtered_result["ra_degrees"] / 15.0
    filtered_result.remove_column("ra_degrees")
    filtered_result["epoch_year"] = EPOCH_YEAR

    result_df = filtered_result.to_pandas(index="hip")
    return result_df


def _get_constellation_boundaries_data_frame(vizier: VizierClass) -> DataFrame:
    DATASET = "VI/49/bound_18"

    constellation_boundary_column_lookup = {
        "cst": "constellation_abbreviation",
        "type": "point_type",
        "_RA.icrs": "ra_degrees",
        "_DE.icrs": "de_degrees",
    }

    results = vizier.query_constraints(
        catalog=DATASET,
        **{"type": "O", "_RA.icrs": "!= NULL", "_DE.icrs": "!= NULL"},
    )
    result = next(iter(results), None)
    if result is None:
        raise RuntimeError()

    filtered_result = result[list(constellation_boundary_column_lookup.keys())]
    for old_name, new_name in constellation_boundary_column_lookup.items():
        filtered_result.rename_column(old_name, new_name)
    filtered_result["ra_hours"] = filtered_result["ra_degrees"] / 15.0
    filtered_result.remove_column("ra_degrees")
    filtered_result_df = filtered_result.to_pandas()

    duplicated_first_vertex_df = (
        filtered_result_df.groupby("constellation_abbreviation", sort=False)
        .apply(lambda group: pd.concat([group, group.head(1)]))
        .reset_index(drop=True)
    )
    duplicated_first_vertex_df["order"] = duplicated_first_vertex_df.groupby(
        "constellation_abbreviation"
    ).cumcount()

    return duplicated_first_vertex_df


def download():
    """Automatically downloads, sanitizes, and organizes astronomical observation data
    from VizieR (https://vizier.cds.unistra.fr/) and Stellarium (https://github.com/Stellarium).
    """

    vizier = Vizier()
    vizier.ROW_LIMIT = -1

    stars_df = _get_stars_data_frame(vizier)
    constellation_boundaries_df = _get_constellation_boundaries_data_frame(vizier)
    pass
