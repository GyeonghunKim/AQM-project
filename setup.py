import io
from setuptools import find_packages, setup


# Read in the README for the long description on PyPI
def long_description():
    with io.open("README.md", "r", encoding="utf-8") as f:
        readme = f.read()
    return readme


setup(
    name="worm_prop",
    version="0.1",
    description="AQM project",
    long_description=long_description(),
    # url="https://github.com/GyeonghunKim/",
    author="",
    author_email="",
    license="",
    packages=find_packages(),
    classifiers=["Programming Language :: Python :: 3"],
    zip_safe=False,
)
