# Star Chart

A simple set of CLI tools designed to automatically download and parse astronomical observation data into a GeoPackage for star chart authoring.

## Installation

Install all dependencies with:

```
poetry install
```

## Usage

All commands listed below should be issued from this project's root folder unless otherwise stated.

### Download data

To download astronomical data for a specific observation date/time (e.g. `2023-06-03` or `2023-06-04T02:00:00Z`) and location (e.g. `-21.232989`, `-44.998945`) run the `download` script as follows:

```sh
poetry run python -m star_chart.commands.download -- 2023-06-04T02:00:00Z -21.232989 -44.998945
```

which is a shorthand for:

```sh
poetry run python -m star_chart.commands.download -c ./config.ini -d ./download -o ./output.gpkg -- 2023-06-04T02:00:00Z -21.232989 -44.998945
```

Given a `./config.ini` file, the script will download all required data to the `./download` fodler, creating a GeoPackage at `./output.gpkg` containing three resulting layers: `stars`, `constellations` and `boundaries`.

Feel free to call `poetry run python -m star_chart.commands.download --help` for help.
