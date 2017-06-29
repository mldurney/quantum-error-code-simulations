#include <stdio.h>
#include <array>
#include "lattices.h"

using namespace std;
const int MAX_FILENAME_SIZE = 255;

void writeOutput(string filename, vector<array<double, 2>> magnetizations);

int main(int argc, char* argv[])
{
    string inFilename;
    int n, updates;
    double t, dt;
    char mode;
    vector<Lattice*> lattices;
    vector<array<double, 2>> magnetizations;

    if (argc == 7)
    {
        inFilename = argv[1];
        t = atof(argv[2]);
        dt = atof(argv[3]);
        n = atoi(argv[4]);
        updates = atoi(argv[5]);
        mode = (argv[6])[0];
    }
    else if (argc == 1)
    {
        cout << "Enter Hamiltonian input file: ";
        cin >> inFilename;
        cout << "Enter minimum temperature (K): ";
        cin >> t;
        cout << "Enter change in temperature between lattices (K): ";
        cin >> dt;
        cout << "Enter total number of lattices: ";
        cin >> n;
        cout << "Enter number of updates per lattice: ";
        cin >> updates;
        cout << "Enter update mode (a - all, r - random): ";
        cin >> mode;
    }
    else
    {
        cout << "Usage: " << argv[0] << " filename(string) ";
        cout << "min_temperature(float) change_temperature(float) ";
        cout << "num_lattices(int) updates(int) mode(char)\n\n";
        exit(1);
    }


    ifstream inFile(inFilename);

    if (!inFile)
    {
        cout << "Invalid file name. " << argv[1] << " does not exist!\n\n";
        exit(EXIT_FAILURE);
    }


    char shape = '\0';
    int rows = -1;
    int cols = -1;

    if (isalpha(inFile.peek()))
    {
        string line;
        string field;
        getline(inFile, line);

        stringstream stream(line);
        stringstream fields[3];

        fields[0] << "\0";
        fields[1] << "-1";
        fields[2] << "-1";

        for (int i = 0; i < 3 && getline(stream, field, ','); ++i)
        {
            fields[i].str(field);
        }

        fields[0] >> shape;
        fields[1] >> rows;
        fields[2] >> cols;
    }


    Hamiltonian h = Hamiltonian(Hamiltonian::importHamiltonian(inFile),
            shape, rows, cols);


    Lattice *nextL;

    for (int i = 0; i < n; i++)
    {
        switch (shape)
        {
            case RECTANGLE:
                nextL = new RectangularLattice(h, t + i * dt, mode);
                break;
            case SQUARE:
                nextL = new SquareLattice(h, t + i * dt, mode);
                break;
            case TRIANGLE:
                nextL = new TriangularLattice(h, t + i * dt, mode);
                break;
            default:
                nextL = new Lattice(h, t + i * dt, mode);
        }

        lattices.push_back(nextL);

        array<double, 2> nextM = {t + i * dt, nextL->getMagnetism()};
        magnetizations.push_back(nextM);
    }


    vector<Lattice*>::iterator itL;
    vector<array<double, 2>>::iterator itM;

    int preupdates = 500;
    for (itL = lattices.begin(); itL != lattices.end(); itL++)
    {
        for (int i = 1; i < preupdates; i++)
        {
            (*itL)->updateLattice();
        }
    }

    for (itL = lattices.begin(), itM = magnetizations.begin();
            itL != lattices.end() && itM != magnetizations.end(); itL++, itM++)
    {
        for (int i = 1; i < updates; i++)
        {
            (*itL)->updateLattice();
            (*itM)[1] += (*itL)->getMagnetism();
        }

        (*itM)[1] = abs((*itM)[1]) / updates;
    }


    string outFilename;

    if (inFilename.find('/') == inFilename.rfind('/'))
    {
        string latticeName = inFilename.substr(0, inFilename.rfind('.'));

        char buffer[MAX_FILENAME_SIZE + 1];
        sprintf(buffer, "magnetizations/%s-%.2fT+%.2fdTx%d-%du-%c.csv",
                latticeName.c_str(), t, dt, n, updates, mode);

        outFilename = buffer;
    }
    else
    {
        string old_dir = "hamiltonians";
        string new_dir = "magnetizations";

        outFilename = inFilename;
        outFilename.replace(outFilename.find(old_dir),
                old_dir.length(), new_dir);
    }

    writeOutput(outFilename, magnetizations);


    for (itL = lattices.begin(); itL != lattices.end(); itL++)
    {
        delete *itL;
    }

    lattices.clear();
    magnetizations.clear();
}

void writeOutput(string filename, vector<array<double, 2>> magnetizations)
{
    bool isNewFile = (ifstream(filename)) ? false : true;
    ofstream file;
    file.open(filename.c_str(), ofstream::out | ofstream::app);

    vector<array<double, 2>>::iterator itM;

    if (isNewFile)
    {
        for (itM = magnetizations.begin(); itM != magnetizations.end(); itM++)
        {
            file << (*itM)[0] << ",";
        }
    }

    file << endl;

    for (itM = magnetizations.begin(); itM != magnetizations.end(); itM++)
    {
        file << (*itM)[1] << ",";
    }

    file.close();

    cout << "Recorded average magnetizations by temperature in " << filename;
    cout << endl;
}
