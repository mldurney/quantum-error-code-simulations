'''
Generate Hamiltonian for lattice of particular shape (rectangle, square,
triangle, or square triangle), given size of lattice and integer coupling that
is same between all neiboring pairs (but may have some probability of being
flipped, lending disorder)
Output Hamiltonian as .csv file
    Format: [shape],rows,cols       <- 's' for square, 'r' for rectangle, etc
            coupling,index1,index2
            1,0,1
            1,0,2
            1,1,2
            ...
'''

import sys
from random import randint
import config as cf


def add_disorder(coupling, disorder):
    '''
    Find coupling when taking disorder into account
    Coupling is original integer coupling of interaction
    Disorder is integer percent chance that coupling will be flipped
    Return opposite coupling at probability given by disorder
    '''

    return coupling if disorder <= randint(0, 99) else coupling * -1

# end add_disorder


def generate_square(neighbors, coupling, disorder, size, coupling2=0):
    '''
    Generate list of interacting indices with their couplings for square
    lattice
    Return Hamiltonian for square lattice of with side length 'size' with same
    coupling between all pairs
    '''

    return generate_rectangle(neighbors, coupling, disorder, size, size, coupling2)

# end generate_square


def generate_rectangle(neighbors, coupling, disorder, rows, cols, coupling2=0):
    '''
    Generate list of interacting indices with their couplings for rectangular
    lattice
    Return Hamiltonian for rectangular lattice of size rows x cols with same
    coupling between all pairs
    '''

    hamiltonian = []

    for index in range(rows * cols):
        # Add interaction between index and right neighbor
        right = (index // cols) * rows + (index + 1) % cols
        hamiltonian.append([add_disorder(coupling, disorder), index, right])

        # Add interaction between index and bottom neighbor
        down = (index + rows) % (rows * cols)
        hamiltonian.append([add_disorder(coupling, disorder), index, down])

        if neighbors == 2:
            # Add interaction between index and second right neighbor
            right2 = (right // cols) * rows + (right + 1) % cols
            hamiltonian.append(
                [add_disorder(coupling2, disorder), index, right2])

            # Add interaction between index and second bottom neighbor
            down2 = (down + rows) % (rows * cols)
            hamiltonian.append(
                [add_disorder(coupling2, disorder), index, down2])

            # Add interaction between index and top-right neighbor
            up_right = (right - rows) % (rows * cols)
            hamiltonian.append(
                [add_disorder(coupling2, disorder), index, up_right])

            # Add interaction between index and bottom-right neighbor
            down_right = (right + rows) % (rows * cols)
            hamiltonian.append(
                [add_disorder(coupling2, disorder), index, down_right])

    return hamiltonian

# end generate_rectangle


def generate_triangle(neighbors, coupling, disorder, rows, colsm, coupling2=0):
    '''
    Generate list of interacting indices with their couplings for triangular
    lattice
    Return Hamiltonian for triangular lattice of size rows x cols with same
    coupling between all pairs
    '''

    hamiltonian = []

    for index in range(rows * cols):
        # Add interaction between index and right neighbor
        right = (index // cols) * rows + (index + 1) % cols
        hamiltonian.append([add_disorder(coupling, disorder), index, right])

        # Add interaction between index and bottom neighbor
        down = (index + rows) % (rows * cols)
        hamiltonian.append([add_disorder(coupling, disorder), index, down])

        # Add interaction between index and bottom-right neighbor
        down_right = (right + rows) % (rows * cols)
        hamiltonian.append(
            [add_disorder(coupling, disorder), index, down_right])

        if neighbors == 2:
            # Add interaction between index and second right neighbor
            right2 = (right // cols) * rows + (right + 1) % cols
            hamiltonian.append(
                [add_disorder(coupling2, disorder), index, right2])

            # Add interaction between index and second bottom neighbor
            down2 = (down + rows) % (rows * cols)
            hamiltonian.append(
                [add_disorder(coupling2, disorder), index, down2])

            # Add interaction between index and 1 down, 2 right neighbor
            down_right2 = (right2 + rows) % (rows * cols)
            hamiltonian.append(
                [add_disorder(coupling2, disorder), index, down_right2])

            # Add interaction between index and top-right neighbor
            up_right = (right - rows) % (rows * cols)
            hamiltonian.append(
                [add_disorder(coupling2, disorder), index, up_right])

            # Add interaction between index and 2 down, 1 right neighbor
            down2_right = (down_right + rows) % (rows * cols)
            hamiltonian.append(
                [add_disorder(coupling2, disorder), index, down_right])

            # Add interactino between index and 2 down, 2 right neighbor
            down2_right2 = (down_right2 + rows) % (rows * cols)
            hamiltonian.append(
                [add_disorder(coupling2, disorder), index, down2_right2])

    return hamiltonian

# end generate_triangle


def generate_striangle(neighbors, coupling, disorder, size, coupling2=0):
    '''
    Generate list of interacting indices with their couplings for square
    triangular lattice
    Return Hamiltonian for square lattice of with side length 'size' with same
    coupling between all pairs
    '''

    return generate_triangle(neighbors, coupling, disorder, size, size, coupling)

# end generate_striangle


def write_hamiltonian(filename, hamiltonian, shape, rows, cols):
    '''
    Receive list of interacting indices with their couplings that make up
    Hamiltonian function for lattice
    Create .csv file ('hamiltonian.csv') giving Hamiltonian function for
    specified lattice
    First line of .csv contains character for shape and row and column count
    Second line contains headers for columns 'coupling,index1,index2'
    Subsequent lines contain each interaction as 'coupling,index1,index2'
    '''

    hamiltonian_file = open(filename, 'w')

    hamiltonian_file.write(shape + ',' + str(rows) + ',' + str(cols) + '\n')
    # hamiltonian_file.write('coupling,index1,index2,rows,cols\n')

    for interaction in hamiltonian:
        hamiltonian_file.write(','.join(map(str, interaction)) + '\n')

    hamiltonian_file.close()

# end write_hamiltonian


def receive_input():
    '''
    Receive user input for lattice shape and size and for coupling between
    each pair of neighboring indices
    Return shape and Hamiltonian function based on input
    '''

    shape = str(input('Enter lattice shape ("s" - square, "r" - rectangle, ' +
                      '"t" - triangle): '))

    neighbors = 2 if int(input('Enter 1 for nearest neighbor interactions, ' +
                               '2 for nearest and next nearest: ')) == 2 else 1

    while shape not in cf.SHAPES:
        shape = str(input('Invalid input. Please enter "s", "r", or "t": '))

    while True:
        try:
            coupling = int(input('Enter integer coupling: '))
            disorder = int(input('Enter disorder [0, 100]: '))
            if neighbors == 2:
                coupling2 = int(
                    input('Enter next-nearest neighbor integer coupling: '))
            else:
                coupling2 = 0
            break
        except ValueError:
            print('Invalid input -- not an integer!')

    if shape in [cf.SQUARE, cf.STRIANGLE]:
        while True:
            try:
                rows = int(input('Enter side length: '))
                cols = rows
                break
            except ValueError:
                print('Invalid input -- not an integer!')

        hamiltonian = (generate_square(neighbors, coupling, disorder, rows, coupling2)
                       if shape == cf.SQUARE
                       else generate_striangle(neighbors, coupling, disorder, rows,
                                               coupling2))

    elif shape in [cf.RECTANGLE, cf.TRIANGLE]:
        while True:
            try:
                rows = int(input('Enter number of rows: '))
                cols = int(input('Enter number of columns: '))
                break
            except ValueError:
                print('Invalid input -- not an integer!')

        hamiltonian = (generate_rectangle(neighbors, coupling, disorder, rows, cols,
                                          coupling2)
                       if shape == cf.RECTANGLE
                       else generate_triangle(neighbors, coupling, disorder, rows, cols,
                                              coupling2))

    return (hamiltonian, shape, rows, cols)
# end receive_input


if __name__ == '__main__':
    FILENAME = sys.argv[1] if len(sys.argv) == 2 else 'hamiltonian'
    write_hamiltonian(FILENAME, *receive_input())
