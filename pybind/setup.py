from setuptools import setup, find_packages

setup(
    name='ACE',
    version='0.0.0',
    install_requires=['numpy'],
    #packages=find_packages(include=['ACE', 'ACE.*']),
    #packages=['ACE'],
    packages=find_packages(),
    package_data={'shared lib': ['*.so']},
    include_package_data=True,
    zip_safe=False,
)
