# extracts file to the current directory
import gzip
import shutil
from typing import IO


def extract_file(downloaded_file: IO):
    compressed_file_name = downloaded_file.name
    extracted_file_name = compressed_file_name.replace(".gz", "")
    with gzip.open(compressed_file_name, "rb") as compressed_file:
        with open(extracted_file_name, "wb") as extracted_file:
            shutil.copyfileobj(compressed_file, extracted_file)
    return extracted_file_name
