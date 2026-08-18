import os
from pathlib import Path
import requests

def download_file(url, output_dir, file_name):
    output_dir = Path(output_dir)

    output_file = output_dir / file_name

    response = requests.get(url, stream=True)
    response.raise_for_status()

    with open(output_file, "wb") as file:
        for chunk in response.iter_content(chunk_size=8192):
            file.write(chunk)

    print(f"Downloaded: {output_file}")

def download_and_extract(url, output_dir, final_filename):
    """Helper function to download and gunzip files."""
    gz_path = os.path.join(output_dir, "temp.gz")
    final_path = os.path.join(output_dir, final_filename)

    if not os.path.exists(final_path):
        print(f"Downloading {final_filename}...")
        response = requests.get(url)
        with open(gz_path, "wb") as f:
            f.write(response.content)
        os.system(f"gunzip -c {gz_path} > {final_path}")
        os.remove(gz_path)

def download_reference_genomes(reference_genomes, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)

    for organism, files in reference_genomes.items():
        for file_type, url in files.items():

            output_name = f"{organism}.{file_type}"

            download_and_extract(
                url,
                output_dir,
                output_name,
            )