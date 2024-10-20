# Star Chart

A simple set of CLI tools designed to:

- Automatically download and parse astronomical observation data into a GeoPackage for star chart authoring.
- Plot star charts using QGIS layouts.

## Installation

Install all dependencies with:

```
uv install
```

## Usage

All commands listed below should be issued from this project's root folder unless otherwise stated.

### `download`

To download astronomical data for a specific observation date/time (e.g. `2023-06-03` or `2023-06-04T02:00:00`) and location (e.g. `-21.232989`, `-44.998945`) run the `download` command as follows:

```sh
python -m star_chart download -- 2023-06-04T02:00:00 -21.232989 -44.998945
```

which is a shorthand for:

```sh
python -m star_chart download -c ./config.ini -d ./download -o ./output.gpkg -- 2023-06-04T02:00:00 -21.232989 -44.998945
```

Given a `./config.ini` file, the script will download all required temporary data to the `./download` fodler, creating a GeoPackage at `./output.gpkg` containing three resulting layers: `stars`, `constellations` and `boundaries`.

For convenience, a sample `config.ini` is available at the root of this project.

Feel free to call `python -m star_chart download --help` for help.

### `plot`

Feel free to call `python -m star_chart plot --help` for help.