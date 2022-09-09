"""Build contact map graphics."""

import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def get_cli_arguments():
    """Parse command line.

    Returns
    -------
    argparse.Namespace
        Object containing arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--contacts",
        action="store",
        type=str,
        help="Input contact file",
        required=True,
    )
    parser.add_argument(
        "--map",
        action="store",
        type=str,
        help="Output contact map image",
        required=True,
    )
    return parser.parse_args()


def create_contact_map(contacts_file, map_file):
    """Create contact map image.

    Parameters
    ----------
    contacts_file : str
        Contact file
    map_file : str
        Contact map image
    """
    contacts = np.loadtxt(contacts_file)
    plt.imshow(np.log2(contacts+1), cmap="hot")
    plt.colorbar()
    plt.savefig(map_file, dpi=300)


if __name__ == "__main__":
    ARGS = get_cli_arguments()
    create_contact_map(ARGS.contacts, ARGS.map)
