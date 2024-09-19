"""Create a HiC-Pro config file from a template."""

import argparse
from pathlib import Path

import jinja2


def is_file(parser, file_path):
    """Check file exists.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Command line argument parser.
    file_path : str
        File path
    Returns
    -------
    str
        File path.
    """
    if not Path(file_path).is_file():
        parser.error(f"The file {file_path} does not exist")
    else:
        return file_path


def get_cli_arguments():
    """Command line argument parser.

    Returns
    -------
    argparse.Namespace
        Object containing arguments.
    """
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument(
        "--template",
        action="store",
        type=lambda name: is_file(parser, name),
        help="Template file for HiC-Pro configuration.",
        required=True,
    )
    required.add_argument(
        "--chromosome-sizes",
        action="store",
        type=lambda name: is_file(parser, name),
        help="File with chromosome sizes.",
        required=True,
    )
    required.add_argument(
        "--genome-fragment",
        action="store",
        type=lambda name: is_file(parser, name),
        help="File with genome restriction sites.",
        required=True,
    )
    required.add_argument(
        "--ligation-site",
        action="store",
        type=str,
        help="Ligation site.",
        required=True,
    )
    required.add_argument(
        "--genome-index-path",
        action="store",
        type=str,
        help="Path to bowtie2 indexes.",
        required=True,
    )
    required.add_argument(
        "--resolutions",
        action="store",
        type=int,
        nargs="+",
        help="HiC resolutions.",
        required=True,
    )
    required.add_argument(
        "--output",
        action="store",
        type=str,
        help="HiC-Pro configuration file.",
        required=True,
    )
    # Add help.
    optional.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    ARGS = get_cli_arguments()

    with open(ARGS.template, "r") as template_file, \
        open(ARGS.output, "w") as out_file:
        template = jinja2.Template(template_file.read())
        out_file.write(
            template.render(
                genome_index_path=str(Path(ARGS.genome_index_path).resolve()) + "/",
                chromosome_sizes=Path(ARGS.chromosome_sizes).resolve(),
                genome_fragment=Path(ARGS.genome_fragment).resolve(),
                ligation_site=ARGS.ligation_site,
                resolutions=" ".join([str(res) for res in ARGS.resolutions])
            )
        )
