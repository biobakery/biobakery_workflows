#!/usr/bin/env python3

import os
import sys
import tarfile
import shutil
import gzip
from pathlib import Path


def get_size(path):
    """Return total size in bytes of a file or all files under a directory."""
    path = Path(path)

    if path.is_file():
        try:
            return path.stat().st_size
        except FileNotFoundError:
            return 0

    total = 0
    for root, dirs, files in os.walk(path):
        for f in files:
            fp = os.path.join(root, f)
            try:
                total += os.path.getsize(fp)
            except FileNotFoundError:
                continue
    return total


def split_directory(input_dir, max_size_bytes):
    """
    Split the contents of input_dir into subfolders named:
    original_dir_chunk_XXX
    """
    input_dir = Path(input_dir)
    chunk_prefix = input_dir.name + "_chunk"

    temp_base = input_dir.parent / f"{input_dir.name}_split_temp"
    shutil.rmtree(temp_base, ignore_errors=True)
    temp_base.mkdir()

    split_index = 1
    current_size = 0
    current_batch = temp_base / f"{chunk_prefix}_{split_index:03d}"
    current_batch.mkdir()

    for item in sorted(input_dir.iterdir()):
        item_size = get_size(item)

        if current_size + item_size > max_size_bytes and current_size > 0:
            split_index += 1
            current_batch = temp_base / f"{chunk_prefix}_{split_index:03d}"
            current_batch.mkdir()
            current_size = 0

        target = current_batch / item.name
        if item.is_dir():
            shutil.copytree(item, target)
        else:
            shutil.copy2(item, target)

        current_size += item_size

    return temp_base, chunk_prefix


def is_already_compressed(path):
    """Return True if file already appears compressed."""
    compressed_suffixes = {
        ".gz", ".bz2", ".xz", ".zip", ".7z", ".rar", ".zst"
    }

    name = path.name.lower()

    if name.endswith(".tar.gz") or name.endswith(".tgz"):
        return True

    return path.suffix.lower() in compressed_suffixes


def gzip_file_in_place(file_path):
    """Gzip file in place and remove original."""
    gz_path = file_path.with_name(file_path.name + ".gz")

    with open(file_path, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    file_path.unlink()


def gzip_all_files(directory):
    """Recursively gzip files except already-compressed ones."""
    for path in directory.rglob("*"):
        if path.is_file() and not is_already_compressed(path):
            print(f"Gzipping {path}")
            gzip_file_in_place(path)


def create_tarballs(split_dir, output_dir):
    """Gzip files inside each chunk directory, then create .tar archives."""
    for part_dir in sorted(split_dir.iterdir()):
        if part_dir.is_dir():
            print(f"Gzipping files inside {part_dir}")
            gzip_all_files(part_dir)

            tar_path = output_dir / f"{part_dir.name}.tar"
            print(f"Archiving {part_dir} → {tar_path}")

            with tarfile.open(tar_path, "w") as tar:
                tar.add(part_dir, arcname=part_dir.name)


def main():
    if len(sys.argv) != 2:
        print("Usage: python split_then_tar.py <input_directory>")
        sys.exit(1)

    input_path = Path(sys.argv[1]).resolve()

    if not input_path.is_dir():
        print(f"Invalid directory: {input_path}")
        sys.exit(1)

    output_dir = input_path.parent / f"{input_path.name}_tars"
    output_dir.mkdir(exist_ok=True)

    max_size_gb = 500
    max_size_bytes = max_size_gb * 1024 ** 3

    print(f"Splitting '{input_path}' into ~{max_size_gb}GB chunks...")
    split_dir, chunk_prefix = split_directory(input_path, max_size_bytes)

    print("Creating .tar archives with individually gzipped files...")
    create_tarballs(split_dir, output_dir)

    print(f"Done! Archives saved to: {output_dir}")

    print("Cleaning up temporary split directories...")
    shutil.rmtree(split_dir)


if __name__ == "__main__":
    main()