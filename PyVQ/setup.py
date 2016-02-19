from distutils.core import setup

setup(
    name='pyvq',
    version='0.1',
    description='Python tools for Virtual Quake/California analysis. Includes Beta versions (not ready for prime-time).',
    author=[
        'mark yoder',
        'kasey schultz'],
    author_email=[
        'mryoder@ucdavis.edu',
        'kwschultz@ucdavis.edu'],
    license="Open Source",
    packages=[
        'pyvq',
        'pyvq.betas'],
    classifiers=["Development Status :: alpha and pre-alpha"])
