#include <bits/stdc++.h>

using namespace std;

typedef unsigned int uint;
typedef long long ll;

#define deg(i) neighbors[i].size()

const std::string WHITESPACE = " \n\r\t\f\v";

std::string ltrim(const std::string &s)
{
    size_t start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}

std::string rtrim(const std::string &s)
{
    size_t end = s.find_last_not_of(WHITESPACE);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

std::string trim(const std::string &s)
{
    return rtrim(ltrim(s));
}

std::vector<int> split(const std::string &s)
{
    std::vector<int> result;
    int j = 0;
    for (int i = 0; i < s.length(); i++)
    {
        if (s[i] == ' ' || s[i] == '\t' || s[i] == '\n' || s[i] == '\r' || s[i] == '\f' || s[i] == '\v')
        {
            result.push_back(stoi(s.substr(j, i - j)));
            j = i;
        }
    }
    if (j < trim(s).length())
        result.push_back(stoi(s.substr(j, s.length() - j)));
    // result.push_back(stoi(s.substr(j, s.length()-j-1)));
    return result;
}

std::string merge(const std::vector<int> &v)
{
    string s = "";
    for (int i : v)
    {
        s += to_string(i) + " ";
    }
    return s;
}

int main(int argc, char *argv[])
{
    assert(argc == 3);
    string graphpath = argv[1];
    string testpath = argv[2];

    vector<vector<int>> lines1;
    vector<vector<int>> lines2;
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
            vector<int> v1 = split(line1);
            vector<int> v2 = split(line2);
            sort(v1.begin(), v1.end());
            sort(v2.begin(), v2.end());
            lines1.push_back(v1);
            lines2.push_back(v2);
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

        sort(lines1.begin(), lines1.end());
        sort(lines2.begin(), lines2.end());
        for (int i = 0; i < lines1.size(); i++)
        {
            if (merge(lines1[i]) != merge(lines2[i]))
            {
                cerr << "Warning: " << merge(lines1[i]) << " != " << merge(lines2[i]) << "\n";
                return 1;
            }
        }

        cout << "SUCCESS: Files verified!\n";
    }
    else
    {
        cerr << "Warning: Unable to open files\n";
        return 1;
    }
    return 0;
}