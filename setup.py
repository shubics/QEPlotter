from pathlib import Path

setup(
    name='qeplotter',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'matplotlib',
        'numpy',
        'seaborn',
        'scipy',
        'pandas',
    ],
    author='Şuayb Yıldız',
    description='Quantum ESPRESSO band structure and DOS plotting tool',
    long_description=Path("README.md").read_text(),
    long_description_content_type='text/markdown',
    url='https://github.com/shubics/QEPlotter',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
)
