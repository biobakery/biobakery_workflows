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

VERSION = "0.2.0"

setuptools.setup(
    name="biobakery_workflows",
    version=VERSION,
    license="MIT",
    description="bioBakery Workflows: A collection of bioBakery AnADAMA2 workflows",
    platforms=['Linux','MacOS'],
    classifiers=[
        "Programming Language :: Python",
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Operating System :: MacOS",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
    install_requires=['anadama2'],
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
