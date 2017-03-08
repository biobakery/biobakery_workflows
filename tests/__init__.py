
import os
import unittest

# find and return all tests in suite
def get_test_suite():
    test_directory=os.path.dirname(os.path.abspath(__file__))
    return unittest.TestLoader().discover(test_directory,pattern='*test*.py')

