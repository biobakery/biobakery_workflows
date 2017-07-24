
import os

def get_file(file):
    # get the full path to the data file
    install_folder=os.path.dirname(os.path.realpath(__file__))
    return os.path.join(install_folder,file)

def get_tutorial_folder():
    return os.path.join(os.path.dirname(os.path.realpath(__file__)),
        os.pardir,os.pardir,"tutorial")

def get_kneaddata_hg_demo_folder():
    return os.path.join(get_tutorial_folder(),"kneaddata_demo_db","Homo_sapiens_demo")

def get_kneaddata_silva_demo_folder():
    return os.path.join(get_tutorial_folder(),"kneaddata_demo_db","SILVA_demo")