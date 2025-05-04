from setuptools import find_packages, setup

setup(
    name='genbenchQC',
    version='0.1',
    packages=find_packages(),
    package_dir={"": "src"},
    tests_require=["pytest"],
    entry_points='''
      [console_scripts]
      evaluate_sequences=genbenchQC.evaluate_sequences:main
      evaluate_dataset=genbenchQC.evaluate_dataset:main
      ''',
)