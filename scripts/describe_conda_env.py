"""Extract conda environment information."""

from watermark import watermark


print(
    watermark(
        machine=True,
        python=True,
        packages="numpy,pandas,scipy,biopandas,matplotlib,Bio,jinja2,sklearn,iced,pastis",
        watermark=True,
        conda=True
    )
)
