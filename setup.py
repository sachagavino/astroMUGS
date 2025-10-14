import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="chemdiskpy",
    version="1.0.1",
    author="Sacha Gavino",
    author_email="sacha.gavino@nbi.ku.dk",
    description="Thermal and chemical modeling of multiple grain-sized protoplanetary disks",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
      packages=[\
        "chemdiskpy",\
        "chemdiskpy.constants", \
        "chemdiskpy.dust",\
        "chemdiskpy.modeling",\
        "chemdiskpy.radmc3d",\
        "chemdiskpy.plotting",\
        "chemdiskpy.nautilus"],\
        package_dir={\
        "chemdiskpy.dust": 'chemdiskpy/dust'}, \
        package_data={\
        'chemdiskpy.nautilus': ['network/*.in']}, \
    python_requires=">=3.6",
)
