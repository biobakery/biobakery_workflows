
import os

def get_file(file):
    # get the full path to the data file
    install_folder=os.path.dirname(os.path.realpath(__file__))
    return os.path.join(install_folder,file)
