import sys
import os


def findDir(dir_name):

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

        for i in range(folders.index(dir_name) + 1):
            dir_path = os.path.join(dir_path, folders[i])

        return dir_path
    
    else:
        raise ValueError('Folder is not a parent directory name!')

# end findDir


SHAPES = [RECTANGLE, SQUARE, TRIANGLE, STRIANGLE] = ['r', 's', 't', 'v']
DIR_NAMES = ['magnetizations', 'binder_cumulants']
ROOT_NAME = 'IsingSimulation'

try:
    ROOT_DIR = findDir(ROOT_NAME)
except ValueError as msg:
    print(msg)
    sys.exit(1)