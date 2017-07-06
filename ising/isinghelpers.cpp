#include "isinghelpers.h"

void openHamiltonianFile(ifstream& file, string filename)
{
    file.open(filename);

    if (!file)
    {
        cout << "Invalid file name. " << filename << " does not exist!\n\n";
        exit(EXIT_FAILURE);
    }
}

Hamiltonian readHamiltonian(ifstream& file, char& shape)
{
    shape = '\0';
    int rows = -1;
    int cols = -1;

    if (isalpha(file.peek()))
    {
        string line;
        string field;
        getline(file, line);

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

    return Hamiltonian(Hamiltonian::importHamiltonianVector(file), shape,
            rows, cols);
}

string getOutFilename(string inFilename, string oldDir, string newDir)
{
    string outFilename;

    if (inFilename.find('/') == inFilename.rfind('/'))
    {
        string objectName = inFilename.substr(0, inFilename.rfind('.'));

        char buffer[MAX_FILENAME_SIZE + 1];
        sprintf(buffer, "%s/%s.csv", newDir.c_str(), objectName.c_str());

        outFilename = buffer;
    }
    else
    {
        outFilename = inFilename;
        outFilename.replace(outFilename.find(oldDir), oldDir.length(), newDir);
    }

    return outFilename;
}

void writeOutput(string filename, vector<double> temp, vector<double> results)
{
    bool isNewFile = (ifstream(filename)) ? false : true;
    ofstream file(filename.c_str(), ofstream::out | ofstream::app);

    vector<double>::iterator it;

    if (isNewFile)
    {
        for (it = temp.begin(); it != temp.end(); ++it)
        {
            file << *it;

            if (it + 1 != temp.end())
            {
                file << ",";
            }
        }

    }

    file << endl;

    for (it = results.begin(); it != results.end(); ++it)
    {
        file << *it;

        if (it + 1 != results.end())
        {
            file << ",";
        }
    }

    file.close();

    cout << "Recorded results by temperature in " << filename;
    cout << endl;
}
