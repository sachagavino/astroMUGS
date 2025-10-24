import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="astromugs",
    version="1.0.5",
    author="Sacha Gavino",
    author_email="sacha.gavino@unibo.it",
    description="Thermal and chemical modeling of multiple grain-sized protoplanetary disks",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
      packages=[\
        "astromugs",\
        "astromugs.constants", \
        "astromugs.dust",\
        "astromugs.modeling",\
        "astromugs.radmc3d",\
        "astromugs.plotting",\
        "astromugs.nautilus"],\
        package_dir={\
        "astromugs.dust": 'astromugs/dust'}, \
        package_data={\
        'astromugs.nautilus': ['network/*.in']}, \
    python_requires=">=3.10",
)
