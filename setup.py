from pathlib import Path

from setuptools import setup, find_packages

ROOT = Path(__file__).resolve().parent
version_namespace = {}
exec(
    (ROOT / "qeplotter" / "version.py").read_text(encoding="utf-8"),
    version_namespace,
)

setup(
    name='qeplotter',
    version=version_namespace["__version__"],
    packages=find_packages(),
    py_modules=['gui_mod'],
    python_requires='>=3.9',
    install_requires=[
        'numpy>=1.24,<3',
        'matplotlib>=3.7,<4',
        'seaborn>=0.12,<1',
        'scipy>=1.10,<2',
        'pandas>=2.0,<3',
        'streamlit==1.50.0',
        'ase>=3.23,<4',
        'spglib>=2.3,<3',
        'plotly>=5.20,<7',
    ],
    author='Şuayb Yıldız',
    author_email='suaybyildiz1@gmail.com',
    description='Quantum ESPRESSO band structure and DOS plotting tool',
    long_description=(ROOT / "README.md").read_text(encoding="utf-8"),
    long_description_content_type='text/markdown',
    url='https://github.com/shubics/QEPlotter',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Operating System :: OS Independent',
    ],
)
