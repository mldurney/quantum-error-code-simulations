import os.path
import pandas as pd
import generate_hamiltonian as gh

class Hamiltonian:
    '''
    TODO: docstring
    '''

    def __init__(self, hamiltonian, shape=None):
        self.hamiltonian = hamiltonian
        self.shape = shape
        self.indices = hamiltonian[index1 : index2].tolist()
        self.local_terms = generate_local_terms(self)

    # end __init__

    def generate_local_terms(self):
        pass
# end Hamiltonian


def main():
    '''
    TODO: docstring
    '''

    if not os.path.isfile(gh.FILENAME):
        gh.write_hamiltonian(*gh.receive_input())

    hamiltonian = pd.read_csv(gh.FILENAME, skiprows=[0])
    shape = pd.read_csv(gh.FILENAME, header=None, nrows=1).iloc[0][0]

    print((hamiltonian['index1']+hamiltonian['index2']).tolist())


if __name__ == '__main__':
    main()
