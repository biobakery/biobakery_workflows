"""
bioBakery Workflows: config module
Configuration settings for workflows and tasks

Copyright (c) 2016 Harvard School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import os
import sys

def get_environ_variable(environ_variable):
    """ Try to get the environment variable """
    
    try:
        variable = os.environ[environ_variable]
    except KeyError:
        sys.exit("ERROR: Please set the environment variable " + environ_variable)
    
    return variable

kneaddata_db_human_genome = get_environ_variable("KNEADDATA_DB_HUMAN_GENOME")
kneaddata_db_human_metatranscriptome = get_environ_variable("KNEADDATA_DB_HUMAN_TRANSCRIPTOME")


