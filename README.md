# Star Chart

A simple set of CLI tools designed to automatically download and parse astronomical observation data into a GeoPackage for star chart authoring.

## Installation

Install all dependencies and scripts with:

```
poetry install
```

## Usage

### Download data

To download astronomical data for a specific observation date (e.g. `2023-06-03`) and location (e.g. `-21.232989`, `-44.998945`) run the `download` script as follows:

```
poetry run python -m star_chart.commands.download -d ./download -o ./output.gpkg -- 2023-06-03 -21.232989 -44.998945 config.ini
```

After downloading all temporary data to the `./download` fodler, a GeoPackage at `./output.gpkg` will be created containing three layers: `stars`, `constellations` and `boundaries`.

Feel free to call `poetry run python -m star_chart.commands.download --help` for help.