#!/usr/bin/env python3
# Set up scratch directory if required.
# Usage: ./scratchSetup.py scratchdir

import os
import sys


def main():
    try:
        if sys.argv[1] != "None":
            os.symlink(sys.argv[1], "Data/Scratch")
        else:
            os.makedirs("Data/Scratch")
    except FileExistsError:
        pass


if __name__ == "__main__":
    main()
