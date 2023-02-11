from setuptools import setup, find_packages

setup(
    name="ao12psr",
    version="0.0.0",
    #packages=find_packages(),
    scripts=["bin/psr12run.py"],
    #package_dir={"ao12psr": "ao12psr"},
    #package_data={"fetch": ["models/model_list.csv", "models/*/*"]},
    #url="https://github.com/devanshkv/fetch",
    #tests_require=["pytest", "pytest-cov"],
    #license="GNU General Public License v3.0",
    author=["Benetge Perera"],
    author_email=["bhakthiperera@gmail.com"],
    description="AO 12 telescope pulsar data processing package",
    classifiers=[
        "Natural Language :: English",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        #"License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
)
