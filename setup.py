import setuptools
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()
with open(path.join(this_directory, 'requirements.txt'), encoding='utf-8') as f:
    requirements = [line.strip() for line in f]

setuptools.setup(
    name="isotools",
    version="0.0.1",
    author="Matthias Lienhard",
    author_email="lienhard@molgen.mpg.de",
    description="framework for the analysis of long read transcriptome sequencing data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    entry_points={
        'console_scripts':['run_isotools = isotools.run_isotools:main']},
    install_requires=requirements
)
