#!/usr/bin/env python3
# Set up scratch directory if required.
# Usage: ./scratchSetup.py scratchdir

import os
import sys


def main():
    try:
        scratchdir = sys.argv[1]
    except IndexError:
        scratchdir = None

    if scratchdir is not None:
        os.symlink(scratchdir, os.getcwd() + "/Data/Scratch")
    else:
        try:
            os.makedirs(os.getcwd() + "/Data/Scratch")
        except FileExistsError:
            pass


if __name__ == "__main__":
    main()
