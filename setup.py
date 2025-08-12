from setuptools import find_packages, setup

requirements = [
    'numpy>=1.23',
    'pandas>=1.5',
    'matplotlib>=3.6',
    'seaborn>=0.12',
    'biopython>=1.8',
    'scikit-learn>=1.2',
    'cdhit-reader==0.2.0',
    'statsmodels>=0.13',
]

test_requirements = [
    'pytest>=3',
]

setup(
    name='genbenchQC',
    version='0.1',
    description='Genomic Benchmarks QC: Automated Quality Control for Genomic Machine Learning Datasets',
    author="Katarina Gresova",
    author_email='gresova11@gmail.com',
    packages=find_packages(),
    package_dir={"": "src"},
    install_requires=requirements,
    extras_require={
        "develop": test_requirements,
    },
    tests_require=["pytest"],
    test_suite='tests',
    entry_points='''
      [console_scripts]
      evaluate_sequences=genbenchQC.evaluate_sequences:main
      evaluate_dataset=genbenchQC.evaluate_dataset:main
      evaluate_split=genbenchQC.evaluate_split:main
      ''',
    keywords=["genomic benchmarks", "deep learning", "machine learning",
      "computational biology", "bioinformatics", "genomics", "quality control"],
)