"""
bioBakery Workflows setup

To run: python setup.py install
"""

import sys

# required python version
required_python_version_major = 2
required_python_version_minor = 7

# Check the python version
try:
    if (sys.version_info[0] != required_python_version_major or
        sys.version_info[1] < required_python_version_minor):
        sys.exit("CRITICAL ERROR: The python version found (version "+
            str(sys.version_info[0])+"."+str(sys.version_info[1])+") "+
            "does not match the version required (version "+
            str(required_python_version_major)+"."+
            str(required_python_version_minor)+"+)")
except (AttributeError,IndexError):
    sys.exit("CRITICAL ERROR: The python version found (version 1) " +
        "does not match the version required (version "+
        str(required_python_version_major)+"."+
        str(required_python_version_minor)+"+)")

try:
    import setuptools
except ImportError:
    sys.exit("Please install setuptools.")
    
from glob import glob    

VERSION = "0.3.1"

AUTHOR = "bioBakery workflows development team"
AUTHOR_EMAIL = "biobakery-users@googlegroups.com"
MAINTAINER = "Lauren McIver"
MAINTAINER_EMAIL = "lauren.j.mciver@gmail.com"

setuptools.setup(
    name="biobakery_workflows",
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    version=VERSION,
    license="MIT",
    description="bioBakery workflows: A collection of meta'omic analysis workflows",
    long_description="bioBakery workflows is a collection of workflows and tasks for "+\
        "executing common microbial community analyses using standardized, validated "+\
        "tools and parameters. Quality control and statistical summary reports are "+\
        "automatically generated for most data types, which include 16S amplicons, "+\
        "metagenomes, and metatranscriptomes. Workflows are run directly from the "+\
        "command line and tasks can be imported to create your own custom workflows. "+\
        "The workflows and tasks are built with AnADAMA2 which allows for parallel "+\
        "task execution locally and in a grid compute environment.",
    url="http://huttenhower.sph.harvard.edu/biobakery_workflows",
    keywords=['microbial','microbiome','bioinformatics','microbiology','metagenomic','metatranscriptomic','anadama2','humann2','metaphlan2','strainphlan'],
    platforms=['Linux','MacOS'],
    classifiers=[
        "Programming Language :: Python",
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Operating System :: MacOS",
        "Operating System :: Unix",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
    install_requires=['anadama2>=0.3.1'],
    packages=setuptools.find_packages(),
    entry_points={
        'console_scripts': [
            'biobakery_workflows = biobakery_workflows.biobakery_workflows:main',
            'biobakery_workflows_databases = biobakery_workflows.biobakery_workflows_databases:main',
        ]},
    package_data={
        'biobakery_workflows' : [
            'workflows/*.py',
            'document_templates/*.py',
            'document_templates/*.mdw',
            'data/*',
        ]},
    data_files=[
        ("tutorial/input/", glob("examples/tutorial/input/*")),
        ("tutorial/kneaddata_demo_db/Homo_sapiens_demo/", glob("examples/tutorial/kneaddata_demo_db/Homo_sapiens_demo/*")),
        ("tutorial/kneaddata_demo_db/SILVA_demo/", glob("examples/tutorial/kneaddata_demo_db/SILVA_demo/*"))
    ],
    scripts=glob('biobakery_workflows/workflows/*py')+glob('biobakery_workflows/scripts/*py'),
    test_suite='tests.get_test_suite',
    zip_safe = False
 )
