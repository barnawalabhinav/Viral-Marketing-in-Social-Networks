#include <bits/stdc++.h>

using namespace std;

typedef unsigned int uint;
typedef long long ll;

#define deg(i) neighbors[i].size()

int main(int argc, char *argv[])
{
    assert(argc == 3);
    string graphpath = argv[1];
    string testpath = argv[2];

    vector<vector<int>> lines;
    ifstream file1(testpath);
    ifstream file2(graphpath);
    if (file1.is_open() && file2.is_open())
    {
        string line1, line2;
        int i = 0;
        while (getline(file1, line1) && getline(file2, line2))
        {
            i++;
            // cout << line1 << "\n";
            // cout << line2 << "\n";
            if (line1 != line2)
            {
                cerr << "Warning: " << line1 << " != " << line2 << "\n";
                file1.close();
                file2.close();
                return 1;
            }
        }
        if (getline(file1, line1) || getline(file2, line2))
        {
            cerr << "Warning: File lengths not equal\n";
            file1.close();
            file2.close();
            return 1;
        }
        file1.close();
        file2.close();

        cout << "SUCCESS: Files verified!\n";
    }
    else
    {
        cerr << "Warning: Unable to open files\n";
        return 1;
    }
    return 0;
}