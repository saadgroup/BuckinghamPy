from setuptools import setup
with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='buckinghampi',
    version='0.1.2',
    packages=[''],
    package_dir={'buckinghampi': 'src'},
    url='',
    license='Apache 2.0',
    author='Mokbel Karam',
    author_email='karammokbel@gmail.com',
    description='Generates the pi terms obtained from Buckingham pi theorem',
    install_requires=required
)
