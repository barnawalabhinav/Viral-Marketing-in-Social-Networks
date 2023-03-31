#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char *argv[])
{
    ifstream inFile;

    if (argc < 4)
    {
        cout << "Usage: " << argv[0] << " <input filename> <output filename> <header filename>" << endl;
        return 1;
    }

    inFile.open(argv[1]);
    if (!inFile)
    {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
    
    ofstream hclean(argv[3], ios::out | ios::binary | ios::trunc);
    hclean.close();
    ofstream hout(argv[3], ios::out | ios::binary | ios::app);

    ofstream fclean(argv[2], ios::out | ios::binary | ios::trunc);
    fclean.close();
    ofstream fout(argv[2], ios::out | ios::binary | ios::app);
    int x, i = 0, j = 2;
    while (inFile >> x)
    {
        fout.write((char*)&x, sizeof(x));
        if (i == j) {
            int k = 4*j;
            hout.write((char*)&k, sizeof(k));
        } else if (i == j + 1) {
            j += x+2;
        }
        i++;
    }

    fout.close();
    hout.close();
    inFile.close();
}