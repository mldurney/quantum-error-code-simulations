#%%
'''
Define common constants (i.e., directory names) for scripts module
'''

import sys
import os


def find_parent_dir(dir_name):
    '''
    Find full directory path of closest parent directory that has basename
    matching dir_name
    Return path (dir_path) if found, otherwise raise ValueError
    '''

    folders = []
    path = os.getcwd()

    while True:
        path, folder = os.path.split(path)

        if folder != '':
            folders.insert(0, folder)

        else:
            if path != '':
                folders.insert(0, path)
            break

    if dir_name in folders:
        dir_path = ''

        for i in range(len(folders) - folders[::-1].index(dir_name)):
            dir_path = os.path.join(dir_path, folders[i])

        return dir_path

    else:
        raise ValueError(dir_name + ' is not a parent directory name!')

# end findDir


def find_exec_path(directory, filename):
    '''
    Find full path for executable file (no extension or .exe) given by
    filename within directory given (may be found in subdirectory)
    Return path (file_path) if found, otherwise raise ValueError
    '''

    valid_filenames = [filename, filename + '.exe']

    for root, _, names in os.walk(directory):
        for name in names:
            if name in valid_filenames:
                return os.path.join(root, name)

    raise ValueError(filename + ' is not an executable file!')

# end find_exec_path


SHAPES = [RECTANGLE, SQUARE, TRIANGLE, STRIANGLE] = ['r', 's', 't', 'v']
DIR_NAMES = ['magnetizations', 'binder_cumulants']
ROOT_NAME = 'IsingSimulation'
MAKE_NAME = 'Makefile'
CSIM_NAME = 'isingsimulation'

try:
    ROOT_DIR = find_parent_dir(ROOT_NAME)
    MAKE_DIR = os.path.dirname(find_exec_path(ROOT_DIR, MAKE_NAME))
except ValueError as msg:
    print('ValueError: ' + str(msg))
    sys.exit(1)
