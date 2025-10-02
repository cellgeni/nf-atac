#!/usr/bin/env python3

import os
import pickle
import glob
import argparse
from pycisTopic.cistopic_class import merge, CistopicObject


def init_parser() -> argparse.ArgumentParser:
    """
    Initialise argument parser for the script
    """
    parser = argparse.ArgumentParser(
        description="Creates cisTopic object from fragments, consensus peaks and quality control files"
    )
    parser.add_argument(
        "--cistopic",
        metavar="<dir>",
        type=str,
        nargs="+",
        help="Specify a path to the directory with cistopic objects",
    )
    parser.add_argument(
        "--output",
        metavar="<file>",
        type=str,
        help="Specify a name to the file with combined cistopic object",
    )
    return parser


def read_cistopic_object(file: str) -> CistopicObject:
    """
    Read cistopic object from the file
    """
    with open(file, "rb") as f:
        cistopic_obj = pickle.load(f)
    return cistopic_obj


def main():
    """
    Main function of the script
    """
    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()

    # combine cistopic objects
    cistopic_objects = [read_cistopic_object(file) for file in args.cistopic]
    cistopic_combined = merge(cistopic_objects)

    # write cistopic object to the file
    with open(args.output, "wb") as f:
        pickle.dump(cistopic_combined, f)


if __name__ == "__main__":
    main()
