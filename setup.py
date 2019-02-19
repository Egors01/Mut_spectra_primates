from setuptools import find_packages, Distribution
from setuptools import setup

def get_module_name(module):
    return module.__name__.split('.')[-1]


class BinaryDistribution(Distribution):
    def is_pure(self):
        return False

def package_files(directory):
    import os
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            if os.path.splitext(filename)[-1] != '.py':
                paths.append(os.path.join(path, filename))
    return paths
extra_files = package_files('primatestools')
# # primary_extra_files = [extra.replace('biota/', '') for extra in primary_extra_files]
#
# package_data.extend(extra_files)
# print(package_data)
#print (find_packages())
setup(
    name='primatestools',
    version='0.3',
    packages=find_packages(),
    url='https://github.com/Egors01/Mut_spectra_primates',
    license='',
    author='egors',
    author_email='none',
    distclass=BinaryDistribution,
    description=''
)
#print(extra_files)
# setup(
#     name='core_biota_secondary',
#     version=0.1,
#     packages=find_packages(),
#     url='',
#     license='',
#     author='Atlas/Knomics',
#     author_email='',
#     description='',
#     include_package_data=True,
#     install_requires=install_requires,
#     distclass=BinaryDistribution,
#     entry_points={
#         'console_scripts': [
#             'run_secondary = biota_secondary.run_secondary:main',
#             'run_biota_module = biota_secondary.modules.run_biota_module:main'
#         ]
#     }
# )