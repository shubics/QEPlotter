from setuptools import setup, find_packages

setup(
    name='qeplotter',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'matplotlib',
        'numpy'
    ],
    author='Şuayb Yıldız',
    description='Quantum ESPRESSO band structure and DOS plotting tool',
    url='https://github.com/shubics/QEPlotter',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
)
