// TODO : Change filesize > Maxfilesize to filesize < Maxfilesize

#include <bits/stdc++.h>
#include <mpi.h>
#include <omp.h>

using namespace std;

typedef unsigned int uint;
typedef long long ll;

#define deg(i) neighbors[i].size()
#define get_node_rank(i) ((i >= (n / size) * (size - 1)) ? size - 1 : i / (n / size))
#define NUM_THREADS 1

int get_edge_rank(int i, int j, int n, int size)
{
    if (i & 1)
    {
        if ((j & 1) && (i > j))
            return get_node_rank(i);
        else
            return get_node_rank(j);
    }
    else
    {
        if (j & 1)
            return get_node_rank(i);
        else if (i < j)
            return get_node_rank(i);
        else
            return get_node_rank(j);
    }
}

int main(int argc, char *argv[])
{
    /***************** PARSE CMD LINE ARGS *****************/

    int taskid = 1;
    int task_p = 1;
    int verbose = 0;
    int startk = -1;
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
        else if (arg.substr(0, 3) == "--p")
            task_p = stoi(arg.substr(4, arg.length() - 4));
        else if (arg == "--help")
        {
            std::cout << "Usage: " << argv[0] << " [--optional value] [--help]\n";
            return 0;
        }
    }
    if (inputpath.length() == 0 || headerpath.length() == 0)
    {
        std::cout << "Optional argument not provided\n";
        abort();
    }
    if (taskid == 2)
        startk = endk;

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided != MPI_THREAD_MULTIPLE)
    {
        printf("MPI_THREAD_MULTIPLE not supported\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /***************** GET INPUT FROM FILE *****************/

    chrono::time_point<std::chrono::system_clock> start, end, start_final, end_final;
    start = std::chrono::system_clock::now();
    start_final = std::chrono::system_clock::now();

    int size_threads = NUM_THREADS;
    int if_present = 0;
    bool proceed = false;
    int totalcount;

    ofstream fout(outputpath, ios::out);

    vector<int> all_deletable;
    vector<int> recvbuf;
    set<pair<int, int>> all_edges;
    vector<set<int>> all_neighbors;
    int *recvcounts = NULL;

#pragma omp parallel num_threads(NUM_THREADS)
    {
        int tid = omp_get_thread_num();

        FILE *ptr = fopen(inputpath.c_str(), "rb"); // r for read, b for binary
        FILE *header = fopen(headerpath.c_str(), "rb");

        unsigned char buffer[4];
        size_t read;

        read = fread(buffer, sizeof(buffer), 1, ptr); // read 8 bytes to our buffer
        int n = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);
        read = fread(buffer, sizeof(buffer), 1, ptr); // read 8 bytes to our buffer
        int m = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);

        set<pair<int, int>> edges;
        vector<set<int>> neighbors(n);
        vector<set<int>> filtered_neighbors(n);

        int degi, node, p, offset, flag, node_bufj, nodej;

        /*************** VERTEX DIVISION OF GRPAH ***************/

        if (rank != size - 1)
        {
            for (int i = rank * (n / size); i < (rank + 1) * (n / size); ++i)
            {
                int index_tid = i - ((rank) * (n / size));

                if (index_tid >= tid * ((n / size) / size_threads) && index_tid < (tid == size_threads - 1 ? (n / size) : ((tid + 1) * ((n / size) / size_threads))))
                {
                    int node = i, p = 0;
                    fseek(header, 4 * node, SEEK_SET);
                    read = fread(buffer, sizeof(buffer), 1, header);

                    offset = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);
                    fseek(ptr, offset + 4, SEEK_SET);
                    read = fread(buffer, sizeof(buffer), 1, ptr);

                    degi = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);

                    unsigned char *buf = (unsigned char *)malloc(degi * sizeof(int));
                    read = fread(buf, degi * sizeof(int), 1, ptr);

                    for (int j = 0; j < degi; ++j)
                    {
                        nodej = (buf[p] & 0xFF) | ((buf[p + 1] & 0xFF) << 8) | ((buf[p + 2] & 0xFF) << 16) | ((buf[p + 3] & 0xFF) << 24);
                        p += 4;
                        (neighbors)[node].insert(nodej);
                        (filtered_neighbors)[node].insert(nodej);

                        // Assign edges to even and smaller (if both even) and larger (if both odd) nodes only
                        if (node & 1)
                        {
                            if ((nodej & 1) && (node > nodej))
                                edges.insert({nodej, node});
                        }
                        else
                        {
                            if (nodej & 1)
                                if (node < nodej)
                                    edges.insert({node, nodej});
                                else
                                    edges.insert({nodej, node});
                            else if (node < nodej)
                                edges.insert({node, nodej});
                        }
                    }
                    std::free(buf);
                }
            }
        }
        else
        {
            int num_verts = n - ((size - 1) * (n / size));
            for (int i = rank * (n / size); i < n; ++i)
            {
                int index_tid = i - ((rank) * (n / size));
                if (index_tid >= tid * ((num_verts) / size_threads) && index_tid < (tid == size_threads - 1 ? (num_verts) : ((tid + 1) * ((num_verts) / size_threads))))
                {
                    int node = i, p = 0;
                    fseek(header, 4 * node, SEEK_SET);
                    read = fread(buffer, sizeof(buffer), 1, header);

                    offset = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);
                    fseek(ptr, offset + 4, SEEK_SET);
                    read = fread(buffer, sizeof(buffer), 1, ptr);

                    degi = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);

                    unsigned char *buf = (unsigned char *)malloc(degi * sizeof(int));
                    read = fread(buf, degi * sizeof(int), 1, ptr);

                    for (int j = 0; j < degi; ++j)
                    {
                        nodej = (buf[p] & 0xFF) | ((buf[p + 1] & 0xFF) << 8) | ((buf[p + 2] & 0xFF) << 16) | ((buf[p + 3] & 0xFF) << 24);
                        p += 4;
                        (neighbors)[node].insert(nodej);
                        (filtered_neighbors)[node].insert(nodej);
                        if (node & 1)
                        {
                            if ((nodej & 1) && (node > nodej))
                                edges.insert({nodej, node});
                        }
                        else
                        {
                            if (nodej & 1)
                                if (node < nodej)
                                    edges.insert({node, nodej});
                                else
                                    edges.insert({nodej, node});
                            else if (node < nodej)
                                edges.insert({node, nodej});
                        }
                    }
                    std::free(buf);
                }
            }
        }

        // /************ FOR INFLUENCER COMPUTATION ************

        /*TODO : add this to influencer in the end*/
        if (taskid == 2 && rank == 0)
        {
        }
        // ************ FOR INFLUENCER COMPUTATION ************/

        // printf("tid = %d, size = %d\n", tid, (int)edges.size());
        for (pair<int, int> edge : edges)
        {
            int node1 = edge.first;
            int node2 = edge.second;

            for (int j = 0; j < 2; ++j)
            {
                int nd = node1;
                if (j > 0)
                    nd = node2;
                if ((neighbors)[nd].empty())
                {
                    fseek(header, 4 * nd, SEEK_SET);
                    read = fread(buffer, sizeof(buffer), 1, header);

                    offset = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);
                    fseek(ptr, offset + 4, SEEK_SET);
                    read = fread(buffer, 4, 1, ptr); // read 8 bytes to our buffer

                    degi = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);

                    unsigned char *buf = (unsigned char *)malloc(degi * sizeof(int));
                    read = fread(buf, degi * sizeof(int), 1, ptr); // read degi * 4 bytes to our buffer

                    p = 0;
                    while (degi--)
                    {
                        nodej = (buf[p] & 0xFF) | ((buf[p + 1] & 0xFF) << 8) | ((buf[p + 2] & 0xFF) << 16) | ((buf[p + 3] & 0xFF) << 24);
                        p += 4;
                        (neighbors)[nd].insert(nodej);
                        (filtered_neighbors)[nd].insert(nodej);
                    }
                    std::free(buf);
                }
            }
        }
        // printf("tid = %d\n", tid);
        std::fclose(ptr);
        std::fclose(header);

        map<pair<int, int>, set<int>> triangles;

        for (pair<int, int> e : edges)
        {
            int i = e.first;
            int j = e.second;
            triangles.insert({{i, j}, set<int>()});
            for (int k : (neighbors)[j])
                if ((neighbors)[i].find(k) != (neighbors)[i].end())
                {
                    triangles.at({i, j}).insert(k);
                    if (i < k)
                        (triangles)[{i, k}].insert(j);
                    else
                        (triangles)[{k, i}].insert(j);

                    if (j < k)
                        (triangles)[{j, k}].insert(i);
                    else
                        (triangles)[{k, j}].insert(i);
                }
        }

        for (int k = startk; k <= endk; k++)
        {
            if (tid == 0)
                if_present = 0;

            while (true)
            {
                auto it = (edges).begin();
                while (it != edges.end())
                {
                    pair<int, int> e = *it;
                    int i = e.first;
                    int j = e.second;
                    if ((int)(triangles.at({i, j})).size() < k)
                    {
#pragma omp critical
                        {
                            all_deletable.push_back(i);
                            all_deletable.push_back(j);
                        }
                        it = edges.erase(it);
                    }
                    else
                        it++;
                }

#pragma omp barrier
                // printf("tid = %d, all_de size = %d, k=%d, edges size = %d\n", tid, all_deletable.size(), k, edges.size());
                if (tid == 0)
                {
                    int sendcount = (int)all_deletable.size(), recvcounts[size], displs[size];
                    MPI_Allgather(&sendcount, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, MPI_COMM_WORLD);
                    proceed = false;
                    for (int i = 0; i < size; i++)
                    {
                        if (recvcounts[i] > 0)
                        {
                            proceed = true;
                            break;
                        }
                    }
                    if (proceed)
                    {

                        displs[0] = 0;
                        for (int i = 1; i < size; i++)
                            displs[i] = displs[i - 1] + recvcounts[i - 1];

                        totalcount = displs[size - 1] + recvcounts[size - 1];

                        recvbuf = vector<int>(totalcount);
                        MPI_Allgatherv(all_deletable.data(), sendcount, MPI_INT, recvbuf.data(), recvcounts, displs, MPI_INT, MPI_COMM_WORLD);

                        all_deletable.clear();
                    }
                }

#pragma omp barrier
                if (proceed == false)
                {
                    break;
                }
                for (int t = 0; t < totalcount; t += 2)
                {
                    int i, j;
                    i = (recvbuf)[t];
                    j = (recvbuf)[t + 1];
                    for (int p : (triangles)[{i, j}])
                    {
                        int w = min(i, p);
                        int x = max(i, p);
                        int y = min(j, p);
                        int z = max(j, p);
                        if ((triangles)[{w, x}].find(j) != (triangles)[{w, x}].end())
                            (triangles)[{w, x}].erase(j);
                        if ((triangles)[{y, z}].find(i) != (triangles)[{y, z}].end())
                            (triangles)[{y, z}].erase(i);
                    }

                    (filtered_neighbors)[i].erase(j);
                    (filtered_neighbors)[j].erase(i);
                }
            }

#pragma omp barrier
            // printf("tid = %d, k = %d\n", tid, k);
            if (taskid == 1 && verbose == 0)
            {
                if (rank == 0)
                {
                    if (((int)(edges).size()) >= 1)
                    {
#pragma omp atomic write
                        if_present = 1;
                    }
#pragma omp barrier
                    // printf("%d, k= %d, rank = %d\n", if_present, k, rank);
                    if (tid == 0)
                    {
                        recvcounts = (int *)malloc(size * sizeof(int));
                        MPI_Gather(&if_present, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
                        int flag = 0;
                        for (int i = 0; i < size; ++i)
                        {
                            // printf("Hi : %d ", recvcounts[i]);
                            if (recvcounts[i] == 1)
                            {
                                fout << 1 << endl;
                                flag = 1;
                                break;
                            }
                        }
                        // printf("\n");
                        if (flag == 0)
                        {
                            fout << 0 << endl;
                        }
                    }
                }
                else
                {
                    if (((int)(edges).size()) >= 1)
                    {
#pragma omp atomic write
                        if_present = 1;
                    }
#pragma omp barrier
                    if (tid == 0)
                        MPI_Gather(&if_present, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
                }
            }
            else
            {
                for (pair<int, int> e : edges)
                {
#pragma omp critical
                    {
                        all_edges.insert(e);
                        cout << e.first << " " << e.second << endl;
                    }
                }
                cout << "----------------- k = " << k << endl;
                for (int i = 0; i < n; i++)
                {
                    for (int j : (filtered_neighbors)[i])
                    {
#pragma omp critical
                        {
                            if (all_neighbors.size() > i)
                                all_neighbors[i].insert(j);
                            else
                            {
                                set<int> tmp;
                                tmp.insert(j);
                                all_neighbors.push_back(tmp);
                            }
                            cout << i << " " << j << "\n";
                        }
                    }
                }
#pragma omp barrier
                printf("HI!\n");
                if (tid == 0)
                {
                    vector<int> head(n, -1);
                    map<int, set<int>> headohead;
                    vector<vector<int>> extr_vert(size);
                    map<int, set<int>> connected_comps;
                    map<int, int> visited;
                    queue<int> trav;
                    auto it = all_edges.begin();
                    set<int> all_heads;

                    while (true)
                    {
                        while (it != all_edges.end() && visited[(*it).first] == 1 && visited[(*it).second] == 1)
                            it++;
                        if (it == all_edges.end())
                            break;

                        int tmp_head = -1;
                        set<int> grp_verts;
                        if (visited[(*it).first] == 0)
                        {
                            trav.push((*it).first);
                            tmp_head = (*it).first;
                            visited[(*it).first] = 1;
                        }
                        if (visited[(*it).second] == 0)
                        {
                            trav.push((*it).second);
                            if (tmp_head == -1)
                                tmp_head = (*it).second;
                            visited[(*it).second] = 1;
                        }
                        if (tmp_head == -1)
                            continue;

                        if (rank == 0)
                            all_heads.insert(tmp_head);

                        while (trav.size() > 0)
                        {
                            int i = trav.front();
                            trav.pop();
                            grp_verts.insert(i);
                            visited[i] = 1;

                            int node_rank = get_node_rank(i);
                            head[i] = tmp_head;
                            if (node_rank != rank)
                            {
                                extr_vert[node_rank].push_back(i);
                                extr_vert[node_rank].push_back(tmp_head);
                            }

                            for (int j : all_neighbors[i])
                            {
                                if (visited[j] == 0)
                                {
                                    trav.push(j);
                                    visited[j] = 1;
                                }
                            }
                        }
                        connected_comps.insert({tmp_head, grp_verts});
                    }

                    int num_stop = rank;
                    vector<int> tmp_conn(2 * n);
                    while (num_stop--)
                    {
                        MPI_Status stat;
                        MPI_Recv(tmp_conn.data(), 2 * n, MPI_INT, MPI_ANY_SOURCE, 5, MPI_COMM_WORLD, &stat);
                        if (tmp_conn[0] != -1)
                        {
                            int cnt = 0;
                            MPI_Get_count(&stat, MPI_INT, &cnt);
                            for (int i = 0; i < cnt; i += 2)
                            {
                                if (head[tmp_conn[i]] == -1)
                                    continue;
                                int hd = head[tmp_conn[i]];
                                headohead[hd].insert(tmp_conn[i + 1]);
                            }
                        }
                    }

                    for (int rk = rank + 1; rk < size; rk++)
                    {
                        if (extr_vert[rk].size() == 0)
                        {
                            int stop = -1;
                            MPI_Send(&stop, 1, MPI_INT, rk, 5, MPI_COMM_WORLD);
                        }
                        else
                        {
                            int sz = extr_vert[rk].size();
                            vector<int> datum = extr_vert[rk];
                            MPI_Send(datum.data(), sz, MPI_INT, rk, 5, MPI_COMM_WORLD);
                        }
                    }

                    if (rank == 0)
                    {
                        int num_stop = size - 1;
                        while (num_stop--)
                        {
                            while (true)
                            {
                                MPI_Status stat;
                                int len;
                                MPI_Recv(&len, 1, MPI_INT, MPI_ANY_SOURCE, 6, MPI_COMM_WORLD, &stat);
                                if (len == -1)
                                    break;
                                int vec[len];
                                int src = stat.MPI_SOURCE;
                                MPI_Recv(vec, len, MPI_INT, src, 7, MPI_COMM_WORLD, &stat);

                                for (int id = 1; id < vec[0]; id++)
                                {
                                    all_heads.insert(vec[id]);
                                    headohead[vec[vec[0]]].insert(vec[id]);
                                    headohead[vec[id]].insert(vec[vec[0]]);
                                }
                                all_heads.insert(vec[1]);
                                for (int ind = vec[0]; ind < len; ind++)
                                {
                                    connected_comps[vec[vec[0]]].insert(vec[ind]);
                                    head[vec[ind]] = vec[vec[0]];
                                }
                            }
                        }

                        set<set<int>> head_comps;
                        map<int, int> visited;
                        queue<int> trav;
                        auto it = all_heads.begin();

                        while (true)
                        {
                            while (it != all_heads.end() && visited[*it] == 1)
                                it++;
                            if (it == all_heads.end())
                                break;

                            set<int> grp_heads;
                            trav.push(*it);
                            visited[*it] = 1;
                            while (trav.size() > 0)
                            {
                                int i = trav.front();
                                trav.pop();
                                grp_heads.insert(i);
                                visited[i] = 1;

                                for (int j : headohead[i])
                                {
                                    if (visited[j] == 0)
                                    {
                                        trav.push(j);
                                        visited[j] = 1;
                                    }
                                }
                            }
                            head_comps.insert(grp_heads);
                        }

                        if (taskid == 1)
                        {
                            map<int, int> vis;
                            if (head_comps.size() == 0)
                                fout << 0 << endl;
                            else
                            {
                                fout << 1 << endl;
                                fout << (int)head_comps.size() << endl;
                                for (auto grp : head_comps)
                                {
                                    for (int hd : grp)
                                    {
                                        for (int v : connected_comps[hd])
                                        {
                                            if (vis[v] != 1)
                                            {
                                                fout << v << " ";
                                                vis[v] = 1;
                                            }
                                        }
                                    }
                                    fout << "\n";
                                }
                            }
                        }
                        else
                        {
                            // /************** CALCULATING INFLUENCER SEQUENTIALLY **************
                            map<int, int> *head_id = new map<int, int>();
                            int id = 0;
                            for (auto grp : head_comps)
                            {
                                for (int hd : grp)
                                    (*head_id)[hd] = id;
                                id++;
                            }

                            vector<set<int>> conn_head_ids(n);
                            vector<int> influencer;
                            for (int vt = 0; vt < n; vt++)
                            {
                                if (head[vt] != -1)
                                    conn_head_ids[vt].insert((*head_id)[head[vt]]);
                                for (int i : neighbors[vt])
                                    if (head[i] != -1)
                                        conn_head_ids[vt].insert((*head_id)[head[i]]);
                                if (conn_head_ids[vt].size() >= task_p)
                                    influencer.push_back(vt);
                            }
                            if (influencer.size() == 0)
                                fout << 0 << endl;
                            else
                            {
                                fout << influencer.size() << endl;
                                int z = head_comps.size();
                                for (int in_vt : influencer)
                                {
                                    if (verbose == 0)
                                        fout << in_vt << " ";
                                    else
                                    {
                                        map<int, int> vis;
                                        fout << in_vt << "\n";
                                        auto itr = head_comps.begin();
                                        for (int i = 0; i < z; i++, itr++)
                                            if (conn_head_ids[in_vt].find(i) != conn_head_ids[in_vt].end())
                                                for (int h : (*itr))
                                                    for (int v : connected_comps[h])
                                                        if (vis[v] != 1)
                                                        {
                                                            fout << v << " ";
                                                            vis[v] = 1;
                                                        }
                                        fout << "\n";
                                    }
                                }
                            }
                            // ************** CALCULATING INFLUENCER SEQUENTIALLY **************/
                        }
                    }
                    else
                    {
                        for (auto comp : connected_comps)
                        {
                            auto it = comp.second.begin();
                            it++;
                            int hdsz = headohead[comp.first].size() + 1;
                            int donesz = 0;
                            bool proceed = true;
                            while (proceed)
                            {
                                int tmpsz = 20000;
                                if (comp.second.size() + hdsz - donesz <= 20000)
                                {
                                    proceed = false;
                                    tmpsz = comp.second.size() + hdsz - donesz;
                                }
                                donesz += tmpsz;
                                vector<int> tmp(tmpsz);
                                int index = 0;
                                tmp[index++] = hdsz;
                                for (auto x : headohead[comp.first])
                                    tmp[index++] = x;
                                tmp[index++] = comp.first;
                                for (; it != comp.second.end() && index < tmpsz; it++)
                                    tmp[index++] = *it;
                                MPI_Send(&index, 1, MPI_INT, 0, 6, MPI_COMM_WORLD);
                                MPI_Send(tmp.data(), index, MPI_INT, 0, 7, MPI_COMM_WORLD);
                            }
                        }
                        int c = -1;
                        MPI_Send(&c, 1, MPI_INT, 0, 6, MPI_COMM_WORLD);
                    }
                    MPI_Barrier(MPI_COMM_WORLD);
                }
            }
#pragma omp barrier
            // printf("tid = %d, here, k = %d, rank = %d\n", tid, k, rank);
        }
        // printf("tid = %d, out of loop\n", tid);
        if (rank == 0)
        {
            if (tid == 0)
            {
                end_final = chrono::system_clock::now();
                chrono::duration<double> elapsed_ms = end_final - start_final;
                std::cout << "Time to end program: " << 1000 * elapsed_ms.count() << " milliseconds\n";
            }
            else
            {
                // cout << "tid = " << tid << '\n';
            }
        }
    }
    MPI_Finalize();
}
