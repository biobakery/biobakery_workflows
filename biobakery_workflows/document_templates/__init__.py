
import os

def get_template(name):
    """ Get the location of the template by name """
    
    # get all of the templates in this folder
    template_extension=".template.py"
    template_install_folder=os.path.dirname(os.path.realpath(__file__))
    templates=filter(lambda file: file.endswith(template_extension),os.listdir(template_install_folder))

    # return the template with the name
    found_template=list(filter(lambda file: file.startswith(name+template_extension), templates))

    if found_template:
        found_template=os.path.join(template_install_folder,found_template[0])
    else:
        found_template=""

    return found_template
