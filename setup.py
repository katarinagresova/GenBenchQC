from setuptools import find_packages, setup

setup(
    name='genData',
    version='0.1',
    packages=find_packages(),
    package_dir={"": "src"},
    tests_require=["pytest"],
    entry_points='''
      [console_scripts]
      evaluate_sequences=genData.evaluate_sequences:main
      evaluate_dataset=genData.evaluate_dataset:main
      generate_negatives=genData.generate_negatives:main
      ''',
)