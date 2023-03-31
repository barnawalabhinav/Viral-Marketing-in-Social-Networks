#include <bits/stdc++.h>

using namespace std;

typedef unsigned int uint;
typedef long long ll;

#define deg(i) neighbors[i].size()

int main(int argc, char *argv[])
{
    /***************** PARSE CMD LINE ARGS *****************/

    int taskid = 1;
    int verbose = 0;
    int startk = 2;
    int endk = 2;
    string inputpath = "";
    string headerpath = "";
    string outputpath = "output.gra";

    for (int i = 1; i < argc; i++)
    {
        string arg = argv[i];
        if (arg.substr(0, 9) == "--verbose")
            verbose = stoi(arg.substr(10, arg.length() - 10));
        else if (arg.substr(0, 8) == "--taskid")
            taskid = stoi(arg.substr(9, arg.length() - 9));
        else if (arg.substr(0, 11) == "--inputpath")
            inputpath = arg.substr(12, arg.length() - 12);
        else if (arg.substr(0, 12) == "--headerpath")
            headerpath = arg.substr(13, arg.length() - 13);
        else if (arg.substr(0, 12) == "--outputpath")
            outputpath = arg.substr(13, arg.length() - 13);
        else if (arg.substr(0, 8) == "--startk")
            startk = stoi(arg.substr(9, arg.length() - 9));
        else if (arg.substr(0, 6) == "--endk")
            endk = stoi(arg.substr(7, arg.length() - 7));
        else if (arg == "--help")
        {
            std::cout << "Usage: " << argv[0] << " [--optional value] [--help]\n";
            return 0;
        }
        else
        {
            std::cerr << "Invalid option: " << arg << "\n";
            return 1;
        }
    }
    if (inputpath.length() == 0 || headerpath.length() == 0)
    {
        std::cout << "Optional argument not provided\n";
        abort();
    }
    if (taskid != 1) // only task 1 is implemented
    {
        std::cout << "Task " << taskid << " not implemented\n";
        abort();
    }

    /***************** GET INPUT FROM FILE *****************/

    chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    unsigned char buffer[8];
    FILE *ptr = fopen(inputpath.c_str(), "rb");          // r for read, b for binary
    size_t read = fread(buffer, sizeof(buffer), 1, ptr); // read 8 bytes to our buffer
    if (read != 1)
    {
        std::cout << "Error reading file\n";
        return 1;
    }
    int n = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);
    int m = (buffer[4] & 0xFF) | ((buffer[5] & 0xFF) << 8) | ((buffer[6] & 0xFF) << 16) | ((buffer[7] & 0xFF) << 24);

    int *graph = (int *)malloc((long long)n * n * sizeof(int)); // adj matrix

    set<pair<int, int>> edges;     // edges
    vector<set<int>> neighbors(n); // adjacency list
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            graph[i * n + j] = -1;

    for (int i = 0; i < n; i++)
    {
        read = fread(buffer, sizeof(buffer), 1, ptr); // read 8 bytes to our buffer
        if (read != 1)
        {
            std::cout << "Error reading file\n";
            return 1;
        }
        int nodei = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);
        int degi = (buffer[4] & 0xFF) | ((buffer[5] & 0xFF) << 8) | ((buffer[6] & 0xFF) << 16) | ((buffer[7] & 0xFF) << 24);
        unsigned char *buf = (unsigned char *)malloc(degi * sizeof(int));

        read = fread(buf, degi * sizeof(int), 1, ptr); // read degi * 4 bytes to our buffer
        if (read != 1)
        {
            std::free(buf);
            std::cout << "Error reading file data" << endl;
            return 1;
        }
        int p = 0;
        for (int i = 0; i < degi; i++)
        {
            int nodej = (buf[p] & 0xFF) | ((buf[p + 1] & 0xFF) << 8) | ((buf[p + 2] & 0xFF) << 16) | ((buf[p + 3] & 0xFF) << 24);
            p += 4;
            graph[nodei * n + nodej] = 0;
            neighbors[nodei].insert(nodej);
            if (nodei < nodej)
                edges.insert({nodei, nodej});
        }
        free(buf);
    }
    std::fclose(ptr);

    end = chrono::system_clock::now();
    chrono::duration<double> elapsed_ms = end - start;
    std::cout << "Time to read input: " << 1000 * elapsed_ms.count() << " milliseconds\n";

    /***************** PRE-PROCESSING GRAPH *****************/

    start = chrono::system_clock::now();

    map<pair<int, int>, set<int>> triangles;
    for (pair<int, int> e : edges)
    {
        int i = e.first;
        int j = e.second;
        for (int k : neighbors[j])
        {
            if (k != i && graph[i * n + k] >= 0)
            {
                graph[i * n + j]++;
                graph[j * n + i]++;
                triangles[{i, j}].insert(k);
            }
        }
    }

    vector<vector<set<int>>> final_ans(endk - startk + 1);
    queue<pair<int, int>> deletable;
    
    for (int k = startk; k <= endk; k++)
    {
        for (pair<int, int> e : edges)
        {
            int i = e.first;
            int j = e.second;
            if (graph[i * n + j] < k)
                deletable.push(e);
        }
        while (deletable.size() > 0)
        {
            pair<int, int> e = deletable.front();
            deletable.pop();
            if (edges.find(e) == edges.end())
                continue;
            int i = e.first;
            int j = e.second;
            edges.erase(e);
            graph[i * n + j] = -1;
            graph[j * n + i] = -1;
            neighbors[i].erase(j);
            neighbors[j].erase(i);
            for (int p : triangles[{i, j}])
            {
                int w = min(i, p);
                int x = max(i, p);
                int y = min(j, p);
                int z = max(j, p);
                triangles[{w, x}].erase(j);
                triangles[{y, z}].erase(i);
                graph[i * n + p]--;
                graph[p * n + i]--;
                graph[j * n + p]--;
                graph[p * n + j]--;
                if (graph[i * n + p] < k)
                    deletable.push({w, x});
                if (graph[j * n + p] < k)
                    deletable.push({y, z});
            }
            triangles.erase(e);
        }
        if (edges.size() == 0)
            continue;

        // for (pair<int, int> e : edges)
        // {
        //     assert(graph[e.first][e.second] >= k);
        //     grp_verts.insert(e.first);
        //     grp_verts.insert(e.second);
        // }

        /****************** K-TRUSS ANALYSIS ******************/

        // for (int i = 0; i < n; i++)
        // {
        //     cout << i << ": ";
        //     for (int j : neighbors[i])
        //         cout << j << " ";
        //     cout << "\n";
        // }
        // for (pair<int, int> e : edges)
        // {
        //     cout << e.first << ", " << e.second << ": " << triangles[e].size() << ": ";
        //     for (int i : triangles[e])
        //         cout << i << " ";
        //     cout << "\n";
        // }

        /****************** K-TRUSS BUILDING ******************/

        map<int, int> visited;
        queue<int> trav;
        auto it = edges.begin();
        while (true)
        {
            while (visited[(*it).first] == 1 && visited[(*it).second] == 1)
                it++;
            if (it == edges.end())
                break;

            set<int> grp_verts;
            if (visited[(*it).first] == 0)
            {
                trav.push((*it).first);
                visited[(*it).first] = 1;
            }
            if (visited[(*it).second] == 0)
            {
                trav.push((*it).second);
                visited[(*it).second] = 1;
            }
            while (trav.size() > 0)
            {
                int i = trav.front();
                trav.pop();
                grp_verts.insert(i);
                visited[i] = 1;
                for (int j : neighbors[i])
                {
                    if (visited[j] == 0)
                    {
                        trav.push(j);
                        visited[j] = 1;
                    }
                }
            }
            final_ans[k - startk].push_back(grp_verts);
        }
    }
    free(graph);

    end = chrono::system_clock::now();
    elapsed_ms = end - start;
    std::cout << "Time to operate: " << 1000 * elapsed_ms.count() << " milliseconds\n";

    /***************** WRITE OUTPUT TO FILE *****************/

    start = chrono::system_clock::now();

    ofstream fclean(outputpath, ios::out | ios::trunc);
    fclean.close();
    ofstream fout(outputpath, ios::out | ios::app);

    for (int i = 0; i <= endk - startk; i++)
    {
        if (final_ans[i].size() == 0)
            fout << 0 << endl;
        else
        {
            fout << 1 << endl;
            if (verbose == 1)
            {
                for (set<int> grp : final_ans[i])
                {
                    int i = 0;
                    for (int v : grp)
                    {
                        fout << v;
                        if (i++ < (int)grp.size() - 1)
                            fout << " ";
                    }
                    fout << endl;
                }
            }
        }
    }
    fout.close();

    end = chrono::system_clock::now();
    elapsed_ms = end - start;
    std::cout << "Time to write output: " << 1000 * elapsed_ms.count() << " milliseconds\n";

    return 0;
}
