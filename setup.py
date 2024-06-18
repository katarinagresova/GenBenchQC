from setuptools import find_packages, setup

setup(
    name='genData',
    version='0.1',
    packages=find_packages(),
    package_dir={"": "src"},
    tests_require=["pytest"],
)