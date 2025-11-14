# astroMUGS

[![astroMUGS](https://img.shields.io/badge/astroMUGS-green?style=flat-square&logo=pypi&logoColor=FFFFFF&labelColor=3A3B3C&color=62F1CD)](https://test.pypi.org/project/astromugs/)
[![image](https://img.shields.io/pypi/pyversions/uv.svg)](https://test.pypi.org/project/astromugs)
[![image](https://img.shields.io/pypi/l/uv.svg)](https://pypi.python.org/pypi/uv)

<picture>
  <source srcset="docs/source/_static/logo-dark.png" media="(prefers-color-scheme: dark)">
  <source srcset="docs/source/_static/logo-light.png" media="(prefers-color-scheme: light)">
  <img src="docs/source/_static/logo-dark.png" alt="astroMUGS", width="300">
</picture>



# A full pipeline for MUlti-Grain Simulations (MUGS) of young stellar objects. 

A modular and intuitive pipeline to couple multi-grain MHD, radiative transfer, and chemistry simulations. 
In the current version, the code converts multi-grain RADMC3D input/output files into a NAUTILUS friendly format, and vice-versa for synthetic line emission maps. 
The pipeline includes ready-to-use sophisticated disk and envelope models to use astroMUGS from scratch. 
In the next release, the pipeline will be based on a combination of [Xarray][1]
and [Zarr][2] for a fast and intuitive use. 
Ultimately, the pipeline will also include post-processing of multi-fluid MHD simulations.

[1]: https://docs.xarray.dev/en/stable/
[2]: https://zarr.dev/


## Installation

- You can find the latest release from test PYPI (for now) [https://test.pypi.org/](https://test.pypi.org/):

        pip install -i https://test.pypi.org/simple/ astromugs


- Alternatively, you can clone the repository:

        git clone https://github.com/sachagavino/astroMUGS.git


## Documentation
A detailed documentation can be found [here][3]. 

[3]: https://astromugs.readthedocs.io/en/latest

## Releases
