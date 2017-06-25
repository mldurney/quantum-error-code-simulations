'''
Generate Hamiltonian for lattice of particular shape (square, rectangle, or
triangle), given size of lattice and integer coupling that is same between
all neiboring pairs
Output Hamiltonian as .csv file ('hamiltonian.csv')
    Format: [shape]                 <- 's' for square, 'r' for rectangle, etc
            coupling,index1,index2
            1,0,1
            1,0,2
            1,1,2
            ...
'''


(SQUARE, RECTANGLE, TRIANGLE) = ('s', 'r', 't')
FILENAME = 'hamiltonian.csv'


def generate_square(coupling, size):
    '''
    Generate list of interacting indices with their couplings for square
    lattice
    Return Hamiltonian for square lattice of with side length 'size' with same
    coupling between all pairs
    '''

    return generate_rectangle(coupling, size, size)

# end generate_square


def generate_rectangle(coupling, rows, cols):
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
        hamiltonian.append([coupling, index, right])

        # Add interaction between index and bottom neighbor
        bottom = (index + rows) % (rows * cols)
        hamiltonian.append([coupling, index, bottom])

    return hamiltonian

# end generate_rectangle


def generate_triangle(coupling, rows, cols):
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
        hamiltonian.append([coupling, index, right])

        # Add interaction between index and bottom neighbor
        bottom = (index + rows) % (rows * cols)
        hamiltonian.append([coupling, index, bottom])

        # Add interaction between index and bottom-right neighbor
        diagonal = (right + rows) % (rows * cols)
        hamiltonian.append([coupling, index, diagonal])

    return hamiltonian

# end generate_triangle


def write_hamiltonian(shape, hamiltonian):
    '''
    Receive list of interacting indices with their couplings that make up
    Hamiltonian function for lattice
    Create .csv file ('hamiltonian.csv') giving Hamiltonian function for
    specified lattice
    First line of .csv contains character for shape
    Second line contains headers for columns 'coupling,index1,index2'
    Subsequent lines contain each interaction as 'coupling,index1,index2'
    '''

    hamiltonian_file = open(FILENAME, 'w')

    hamiltonian_file.write(shape + '\n')
    hamiltonian_file.write('coupling,index1,index2\n')

    for interaction in hamiltonian:
        hamiltonian_file.write(','.join(map(str, interaction)) + '\n')

    hamiltonian_file.close()

#end write_hamiltonian


def receive_input():
    '''
    Receive user input for lattice shape and size and for coupling between
    each pair of neighboring indices
    Return shape and Hamiltonian function based on input
    '''

    shape = str(input('Enter lattice shape ("s" - square, "r" - rectangle, \
                      "t" - triangle): '))

    while (shape not in [SQUARE, RECTANGLE, TRIANGLE]):
        shape = str(input('Invalid input. Please enter "s", "r", or "t": '))


    while True:
        try:
            coupling = int(input('Enter integer coupling: '))
            break
        except ValueError:
            print('Invalid input -- not an integer!')


    if shape == SQUARE:
        while True:
            try:
                size = int(input('Enter square side length: '))
                break
            except ValueError:
                print('Invalid input -- not an integer!')

        hamiltonian = generate_square(coupling, size)

    elif shape == RECTANGLE:
        while True:
            try:
                rows = int(input('Enter number of rows: '))
                cols = int(input('Enter number of columns: '))
                break
            except ValueError:
                print('Invalid input -- not an integer!')

        hamiltonian = generate_rectangle(coupling, rows, cols)

    elif shape == TRIANGLE:
        while True:
            try:
                rows = int(input('Enter number of rows: '))
                cols = int(input('Enter number of columns: '))
                break
            except ValueError:
                print('Invalid input -- not an integer')

        hamiltonian = generate_triangle(coupling, rows, cols)


    return (shape, hamiltonian)

# end main


if __name__ == '__main__':
    write_hamiltonian(*receive_input())