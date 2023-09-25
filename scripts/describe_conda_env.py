"""Extract conda environment information."""

from watermark import watermark


def list_versions():
    """List package versions and describe environment."""
    print(
        watermark(
            machine=True,
            python=True,
            packages="numpy,pandas,scipy,biopandas,matplotlib,Bio,jinja2,sklearn,iced,pastis",
            watermark=True,
            conda=True
        )
    )


if __name__ == "__main__":
    print("-" * 40)
    list_versions()
    print("-" * 40)