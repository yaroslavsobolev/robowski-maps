from robowski.settings import *
import importlib_resources
import os
repo_data_path = str(importlib_resources.files("robowski")) + '/'
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'