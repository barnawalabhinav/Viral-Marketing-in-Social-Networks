// TODO : Change filesize > Maxfilesize to filesize < Maxfilesize

#include <bits/stdc++.h>
#include <mpi.h>
#include <omp.h>

using namespace std;

typedef unsigned int uint;
typedef long long ll;

#define deg(i) neighbors[i].size()
#define get_node_rank(i) ((i >= (n / size) * (size - 1)) ? size - 1 : i / (n / size))

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

    FILE *ptr = fopen(inputpath.c_str(), "rb"); // r for read, b for binary
    FILE *header = fopen(headerpath.c_str(), "rb");

    unsigned char buffer[4];
    size_t read;

    read = fread(buffer, sizeof(buffer), 1, ptr); // read 8 bytes to our buffer
    int n = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);
    read = fread(buffer, sizeof(buffer), 1, ptr); // read 8 bytes to our buffer
    int m = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);

    set<pair<int, int>> *edges = new set<pair<int, int>>();
    vector<set<int>> *neighbors = new vector<set<int>>(n);
    vector<set<int>> *filtered_neighbors = new vector<set<int>>(n);

    int degi, node, p, offset, flag, node_bufj, nodej;

    /*************** VERTEX DIVISION OF GRPAH ***************/

    if (rank != size - 1)
    {
        for (int i = rank * (n / size); i < (rank + 1) * (n / size); ++i)
        {
            int node = i, p = 0;
            // printf("rank = %d has node = %d\n", rank, node);
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
                (*neighbors)[node].insert(nodej);
                (*filtered_neighbors)[node].insert(nodej);
                // Assign edges to even and smaller (if both even) and larger (if both odd) nodes only
                if (node & 1)
                {
                    if ((nodej & 1) && (node > nodej))
                        edges->insert({nodej, node});
                }
                else
                {
                    if (nodej & 1)
                        if (node < nodej)
                            edges->insert({node, nodej});
                        else
                            edges->insert({nodej, node});
                    else if (node < nodej)
                        edges->insert({node, nodej});
                }
                // if (node < nodej)
                //     edges->insert({node, nodej});
            }
            std::free(buf);
        }
    }
    else
    {
        for (int i = rank * (n / size); i < n; ++i)
        {
            int node = i, p = 0;
            // printf("rank = %d has node = %d\n", rank, node);
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
                (*neighbors)[node].insert(nodej);
                (*filtered_neighbors)[node].insert(nodej);
                if (node & 1)
                {
                    if ((nodej & 1) && (node > nodej))
                        edges->insert({nodej, node});
                }
                else
                {
                    if (nodej & 1)
                        if (node < nodej)
                            edges->insert({node, nodej});
                        else
                            edges->insert({nodej, node});
                    else if (node < nodej)
                        edges->insert({node, nodej});
                }
                // if (node < nodej)
                //     edges->insert({node, nodej});
            }
            std::free(buf);
        }
    }

    // if (size == 1)
    // {
    //     for (int i = 0; i < n; ++i)
    //     {
    //         node = i, p = 0;
    //         fseek(header, 4 * node, SEEK_SET);
    //         read = fread(buffer, sizeof(buffer), 1, header);

    //         offset = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);
    //         fseek(ptr, offset + 4, SEEK_SET);
    //         read = fread(buffer, sizeof(buffer), 1, ptr);

    //         degi = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);

    //         unsigned char *buf = (unsigned char *)malloc(degi * sizeof(int));
    //         read = fread(buf, degi * sizeof(int), 1, ptr);

    //         int flag = 0, node_bufj;
    //         for (int j = 0; j < degi; ++j)
    //         {
    //             nodej = (buf[p] & 0xFF) | ((buf[p + 1] & 0xFF) << 8) | ((buf[p + 2] & 0xFF) << 16) | ((buf[p + 3] & 0xFF) << 24);
    //             p += 4;

    //             if (node < nodej)
    //             {
    //                 if (edges->find({node, nodej}) == edges->end())
    //                 {
    //                     edges->insert({node, nodej});
    //                 }
    //             }
    //         }
    //         std::free(buf);
    //     }
    // }
    // else
    // {
    //     if (rank == 0)
    //     {
    //         for (int i = 0; i < n; ++i)
    //         {
    //             node = i, p = 0;
    //             fseek(header, 4 * node, SEEK_SET);
    //             read = fread(buffer, sizeof(buffer), 1, header);

    //             offset = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);
    //             fseek(ptr, offset + 4, SEEK_SET);
    //             read = fread(buffer, sizeof(buffer), 1, ptr);

    //             degi = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);

    //             unsigned char *buf = (unsigned char *)malloc(degi * sizeof(int));
    //             read = fread(buf, degi * sizeof(int), 1, ptr);

    //             int flag = 0, node_bufj;
    //             for (int j = 0; j < degi; ++j)
    //             {
    //                 nodej = (buf[p] & 0xFF) | ((buf[p + 1] & 0xFF) << 8) | ((buf[p + 2] & 0xFF) << 16) | ((buf[p + 3] & 0xFF) << 24);
    //                 p += 4;

    //                 if (node < nodej)
    //                 {
    //                     if (edges->find({node, nodej}) == edges->end())
    //                     {
    //                         edges->insert({node, nodej});

    //                         if ((int)edges->size() == m / size)
    //                         {
    //                             flag = 1;
    //                             node_bufj = p;
    //                             break;
    //                         }
    //                     }
    //                 }
    //             }
    //             std::free(buf);
    //             if (flag == 1)
    //             {
    //                 MPI_Send(&node, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    //                 MPI_Send(&node_bufj, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    //                 break;
    //             }
    //         }
    //     }
    //     else if (rank == size - 1)
    //     {
    //         MPI_Recv(&node, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         MPI_Recv(&node_bufj, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    //         p = node_bufj;
    //         for (int i = node; i < n; ++i)
    //         {
    //             node = i;
    //             fseek(header, 4 * node, SEEK_SET);
    //             read = fread(buffer, sizeof(buffer), 1, header);

    //             offset = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);
    //             fseek(ptr, offset + 4, SEEK_SET);
    //             read = fread(buffer, sizeof(buffer), 1, ptr);

    //             degi = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);

    //             unsigned char *buf = (unsigned char *)malloc(degi * sizeof(int));
    //             read = fread(buf, degi * sizeof(int), 1, ptr);

    //             flag = 0;

    //             for (int j = p / 4; j < degi; ++j)
    //             {
    //                 nodej = (buf[p] & 0xFF) | ((buf[p + 1] & 0xFF) << 8) | ((buf[p + 2] & 0xFF) << 16) | ((buf[p + 3] & 0xFF) << 24);
    //                 p += 4;

    //                 if (node < nodej)
    //                 {
    //                     if (edges->find({node, nodej}) == edges->end())
    //                     {
    //                         edges->insert({node, nodej});
    //                     }
    //                 }
    //             }
    //             std::free(buf);
    //             p = 0;
    //         }
    //     }
    //     else
    //     {
    //         MPI_Recv(&node, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         MPI_Recv(&node_bufj, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    //         p = node_bufj;
    //         for (int i = node; i < n; ++i)
    //         {
    //             node = i;
    //             fseek(header, 4 * node, SEEK_SET);
    //             read = fread(buffer, sizeof(buffer), 1, header);

    //             offset = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);
    //             fseek(ptr, offset + 4, SEEK_SET);
    //             read = fread(buffer, sizeof(buffer), 1, ptr);

    //             degi = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);

    //             unsigned char *buf = (unsigned char *)malloc(degi * sizeof(int));
    //             read = fread(buf, degi * sizeof(int), 1, ptr);

    //             flag = 0;
    //             for (int j = p / 4; j < degi; ++j)
    //             {
    //                 nodej = (buf[p] & 0xFF) | ((buf[p + 1] & 0xFF) << 8) | ((buf[p + 2] & 0xFF) << 16) | ((buf[p + 3] & 0xFF) << 24);
    //                 p += 4;

    //                 if (node < nodej)
    //                 {
    //                     if (edges->find({node, nodej}) == edges->end())
    //                     {
    //                         edges->insert({node, nodej});
    //                         if ((int)edges->size() == m / size)
    //                         {
    //                             flag = 1;
    //                             node_bufj = p;
    //                             break;
    //                         }
    //                     }
    //                 }
    //             }
    //             std::free(buf);
    //             p = 0;
    //             if (flag == 1)
    //             {
    //                 MPI_Send(&node, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    //                 MPI_Send(&node_bufj, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    //                 break;
    //             }
    //         }
    //     }
    // }
    // vector<set<int>> *neighbors = new vector<set<int>>(n);
    // vector<bool> *visit = new vector<bool>(n, false);

    // /************ FOR SEQUENTIAL INFLUENCER COMPUTATION ************
    if (taskid == 2 && rank == 0)
    {
        for (int nd = (n / size); nd < n; ++nd)
        {
            int tmp_p = 0;
            // printf("rank = %d has node = %d\n", rank, node);
            fseek(header, 4 * nd, SEEK_SET);
            read = fread(buffer, sizeof(buffer), 1, header);

            offset = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);
            fseek(ptr, offset + 4, SEEK_SET);
            read = fread(buffer, sizeof(buffer), 1, ptr);

            degi = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);

            unsigned char *buf = (unsigned char *)malloc(degi * sizeof(int));
            read = fread(buf, degi * sizeof(int), 1, ptr);

            for (int j = 0; j < degi; ++j)
            {
                nodej = (buf[tmp_p] & 0xFF) | ((buf[tmp_p + 1] & 0xFF) << 8) | ((buf[tmp_p + 2] & 0xFF) << 16) | ((buf[tmp_p + 3] & 0xFF) << 24);
                tmp_p += 4;
                (*neighbors)[nd].insert(nodej);
            }
            free(buf);
        }
    }
    // ************ FOR SEQUENTIAL INFLUENCER COMPUTATION ************/

    for (pair<int, int> edge : *edges)
    {
        int node1 = edge.first;
        int node2 = edge.second;

        for (int j = 0; j < 2; ++j)
        {
            int nd = node1;
            if (j > 0)
                nd = node2;
            if ((*neighbors)[nd].empty())
            {
                // visit->at(nd) = true;
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
                    (*neighbors)[nd].insert(nodej);
                    (*filtered_neighbors)[nd].insert(nodej);
                }
                std::free(buf);
            }
        }
    }
    std::fclose(ptr);
    std::fclose(header);

    // end = chrono::system_clock::now();
    // chrono::duration<double> elapsed_ms = end - start;
    // std::cout << "Time to read input: " << 1000 * elapsed_ms.count() << " milliseconds\n";

    /***************** PRE-PROCESSING GRAPH *****************/

    // start = chrono::system_clock::now();

    map<pair<int, int>, set<int>> *triangles = new map<pair<int, int>, set<int>>();

    for (pair<int, int> e : *edges)
    {
        int i = e.first;
        int j = e.second;
        triangles->insert({{i, j}, set<int>()});
        for (int k : (*neighbors)[j])
            if ((*neighbors)[i].find(k) != (*neighbors)[i].end())
            {
                triangles->at({i, j}).insert(k);
                if (i < k)
                    (*triangles)[{i, k}].insert(j);
                else
                    (*triangles)[{k, i}].insert(j);

                if (j < k)
                    (*triangles)[{j, k}].insert(i);
                else
                    (*triangles)[{k, j}].insert(i);
            }
    }

    // vector<set<set<int>>> *final_ans = new vector<set<set<int>>>(endk - startk + 1);

    vector<int> *deletable = new vector<int>();

    ofstream fout(outputpath, ios::out);

    // bool finished = false;
    for (int k = startk; k <= endk; k++)
    {
        // cout << k << endl;
        // printf("%d: Edges size = %d\n", k, edges->size());
        // chrono::duration<double> g_time = chrono::system_clock::now() - chrono::system_clock::now();
        while (true)
        {
            // std::printf("rank %d has k = %d\n", rank, k);

            auto it = (*edges).begin();
            while (it != edges->end())
            {
                pair<int, int> e = *it;
                int i = e.first;
                int j = e.second;
                if ((int)(triangles->at({i, j})).size() < k)
                {
                    deletable->push_back(i);
                    deletable->push_back(j);
                    it = edges->erase(it);
                }
                else
                    it++;
            }
            int sendcount = (int)deletable->size(), recvcounts[size], displs[size], totalcount = 0;

            // printf("size %d\n", size);

            // chrono::time_point<std::chrono::system_clock> gather_beg = chrono::system_clock::now();
            //  MPI_Barrier(MPI_COMM_WORLD);
            MPI_Allgather(&sendcount, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, MPI_COMM_WORLD);
            // MPI_Barrier(MPI_COMM_WORLD);

            bool proceed = false;
            for (int i = 0; i < size; i++)
            {
                // if (rank == 0)
                //     printf("%d: %d \n", k, recvcounts[i]);

                if (recvcounts[i] > 0)
                {
                    proceed = true;
                    break;
                }
            }

            // printf("proceed %d\n", proceed);

            if (!proceed)
                break;

            displs[0] = 0;
            for (int i = 1; i < size; i++)
            {
                displs[i] = displs[i - 1] + recvcounts[i - 1];
            }
            totalcount = displs[size - 1] + recvcounts[size - 1];

            // printf("totalcount %d\n", totalcount);

            vector<int> *recvbuf = new vector<int>(totalcount);

            // MPI_Barrier(MPI_COMM_WORLD);
            MPI_Allgatherv(deletable->data(), sendcount, MPI_INT, recvbuf->data(), recvcounts, displs, MPI_INT, MPI_COMM_WORLD);
            // MPI_Barrier(MPI_COMM_WORLD);
            // chrono::time_point<std::chrono::system_clock> gather_end = chrono::system_clock::now();
            // g_time += gather_end - gather_beg;

            deletable->clear();
            // vector<pair<int, int>> *deleted = new vector<pair<int, int>>();
            // for (int i = 0; i < totalcount; i += 2)
            //     deleted->push_back({recvbuf->at(i), recvbuf->at(i + 1)});

            // printf("deleted size %d\n", deleted->size());

            for (int t = 0; t < totalcount; t += 2)
            {
                int i = (*recvbuf)[t];
                int j = (*recvbuf)[t + 1];
                for (int p : (*triangles)[{i, j}])
                {
                    int w = min(i, p);
                    int x = max(i, p);
                    int y = min(j, p);
                    int z = max(j, p);
                    if ((*triangles)[{w, x}].find(j) != (*triangles)[{w, x}].end())
                        (*triangles)[{w, x}].erase(j);
                    if ((*triangles)[{y, z}].find(i) != (*triangles)[{y, z}].end())
                        (*triangles)[{y, z}].erase(i);
                    // graph[i * n + p]--;
                    // graph[p * n + i]--;
                    // graph[j * n + p]--;
                    // graph[p * n + j]--;
                    // if (graph[i * n + p] < k)
                    //     deletable.push({w, x});
                    // if (graph[j * n + p] < k)
                    //     deletable.push({y, z});
                }

                (*filtered_neighbors)[i].erase(j);
                (*filtered_neighbors)[j].erase(i);

                // for (int p : (*neighbors)[i])
                // {
                //     if (p == j)
                //         continue;
                //     int w = min(i, p);
                //     int x = max(i, p);
                //     if ((*neighbors)[p].find(j) != (*neighbors)[p].end())
                //     {
                //         // std::printf("R %d: %d, %d: %d\n", rank, i, j, p);
                //         if ((*triangles)[{w, x}].find(j) != (*triangles)[{w, x}].end())
                //             (*triangles)[{w, x}].erase(j);
                //     }
                // }
                // for (int p : (*neighbors)[j])
                // {
                //     if (p == i)
                //         continue;
                //     int w = min(j, p);
                //     int x = max(j, p);
                //     if ((*neighbors)[p].find(i) != (*neighbors)[p].end())
                //     {
                //         // std::printf("R %d: %d, %d: %d\n", rank, i, j, p);
                //         if ((*triangles)[{w, x}].find(i) != (*triangles)[{w, x}].end())
                //             (*triangles)[{w, x}].erase(i);
                //     }
                // }
            }
        }

        // printf("Time to gather in %d truss: %f milliseconds\n", k, 1000 * g_time.count());

        /************ TEST FOR THE OUTPUT CORRECTNESS BEGINS ************/

        // printf("Edges size = %d\n", edges->size());
        // MPI_Barrier(MPI_COMM_WORLD);

        // set<int> *grp_verts = new set<int>();

        // for (pair<int, int> e : (*edges))
        // {
        //     grp_verts->insert(e.first);
        //     grp_verts->insert(e.second);
        // }
        // // printf("grp_verts size : %d\n", (int)grp_verts->size());

        // // for (int i : *grp_verts)
        // //     printf("%d\n", i);

        // vector<int> group(grp_verts->begin(), grp_verts->end());
        // int sendcount = (int)group.size(), recvcounts[size], displs[size], totalcount = 0;

        // // printf("size %d\n", group.size());

        // MPI_Barrier(MPI_COMM_WORLD);
        // MPI_Gather(&sendcount, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
        // MPI_Barrier(MPI_COMM_WORLD);

        // vector<int> *recvbuf;

        // if (rank == 0)
        // {
        //     displs[0] = 0;
        //     for (int i = 1; i < size; i++)
        //         displs[i] = displs[i - 1] + recvcounts[i - 1];
        //     totalcount = displs[size - 1] + recvcounts[size - 1];
        //     // printf("totalcount %d\n", totalcount);
        //     recvbuf = new vector<int>(totalcount);
        // }

        // MPI_Barrier(MPI_COMM_WORLD);
        // MPI_Gatherv(group.data(), sendcount, MPI_INT, recvbuf->data(), recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
        // MPI_Barrier(MPI_COMM_WORLD);

        // // set<int> *grp_verts_all = new set<int>(recvbuf->begin(), recvbuf->end());
        // // if (grp_verts_all->size() == 0)
        // //     finished = true;

        // if (rank == 0)
        // {
        //     set<int> *grp_verts_all = new set<int>(recvbuf->begin(), recvbuf->end());
        //     if (grp_verts_all->size() == 0)
        //     {
        //         printf("0\n");
        //     }
        //     else
        //     {
        //         printf("1\n");
        //         int j = grp_verts_all->size() - 1;
        //         for (int i : (*grp_verts_all))
        //         {
        //             printf("%d", i);
        //             if (j-- > 0)
        //                 printf(" ");
        //         }
        //         printf("\n");
        //     }
        // }

        /************ TEST FOR THE OUTPUT CORRECTNESS ENDS ************/

        // end = chrono::system_clock::now();
        // elapsed_ms = end - start;
        // std::printf("Time to find %d-truss: %f milliseconds\n", k, 1000 * elapsed_ms.count());

        if (taskid == 1 && verbose == 0)
        {
            int *recvcounts = NULL;
            if (rank == 0)
            {
                int if_present = 0;
                if (((int)(*edges).size()) >= 1)
                    if_present = 1;
                else
                    if_present = 0;

                recvcounts = (int *)malloc(size * sizeof(int));
                MPI_Gather(&if_present, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
                int flag = 0;
                for (int i = 0; i < size; ++i)
                {
                    if (recvcounts[i] == 1)
                    {
                        fout << 1 << endl;
                        flag = 1;
                        break;
                    }
                }
                if (flag == 0)
                {
                    fout << 0 << endl;
                }
            }
            else
            {
                int if_present = 0;
                if (((int)(*edges).size()) >= 1)
                    if_present = 1;
                else
                    if_present = 0;

                MPI_Gather(&if_present, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
            }
        }
        else
        {
            // int last_node = (rank == size - 1) ? n : (rank + 1) * (n / size);
            // int first_node = rank * (n / size);
            // int num_nodes = last_node - first_node;
            // vector<int> *head = new vector<int>(num_nodes, -1);

            ifstream file(inputpath.c_str(), ios::binary);
            file.seekg(0, ios::end);
            int filesize = file.tellg();

            int MAXFILESIZE = 0.5 * (1000000000);

            // if (filesize < MAXFILESIZE)
            if (false)
            {
                // start = chrono::system_clock::now();

                // printf("hello1\n");

                // set<int> *grp_verts = new set<int>();

                // for (pair<int, int> e : (*edges))
                // {
                //     grp_verts->insert(e.first);
                //     grp_verts->insert(e.second);
                // }

                // for(auto v:(*grp_verts))
                // {
                //     printf("rank %d has vert = %d", rank, v);
                // }
                // printf("\n");

                vector<int> *edge_to_vert = new vector<int>();
                for (auto e : *edges)
                {
                    edge_to_vert->push_back(e.first);
                    edge_to_vert->push_back(e.second);
                }
                int sendcount = (int)edge_to_vert->size();
                int *recvcounts = NULL;
                int *displs = NULL;
                int totalcount = 0;

                int *recvbuf = NULL;

                // printf("hello1\n");

                map<int, set<int>> *connected_comps = new map<int, set<int>>();

                if (rank != 0)
                {
                    MPI_Gather(&sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    // printf("rank %d has %d count\n", rank, sendcount);
                    MPI_Gatherv(edge_to_vert->data(), sendcount, MPI_INT, recvbuf, recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
                    // printf("rank %d has %d count\n", rank, sendcount);
                    // if (taskid == 2)
                    // {
                    //     int node_head[2];
                    //     while (true)
                    //     {
                    //         MPI_Recv(&node_head, 2, MPI_INT, 0, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    //         if (node_head[0] == -1)
                    //             break;
                    //         (*head)[node_head[0] - first_node] = node_head[1];
                    //     }
                    // }
                }
                else
                {

                    vector<int> *head = new vector<int>(n, -1);

                    recvcounts = (int *)malloc(size * sizeof(int));

                    MPI_Gather(&sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    // MPI_Barrier(MPI_COMM_WORLD);
                    // printf("rank %d has %d count\n", rank, sendcount);

                    displs = (int *)malloc(size * sizeof(int));

                    displs[0] = 0;

                    for (int i = 1; i < size; i++)
                    {
                        displs[i] = displs[i - 1] + recvcounts[i - 1];
                    }
                    totalcount = displs[size - 1] + recvcounts[size - 1];

                    recvbuf = (int *)malloc(totalcount * sizeof(int));
                    // printf("totalcount %d\n", totalcount);

                    // MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Gatherv(edge_to_vert->data(), sendcount, MPI_INT, recvbuf, recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
                    // printf("rank %d has %d count\n", rank, sendcount);
                    //  MPI_Barrier(MPI_COMM_WORLD);
                    // chrono::time_point<std::chrono::system_clock> gather_end = chrono::system_clock::now();
                    // g_time += gather_end - gather_beg;

                    vector<pair<int, int>> *all_edges = new vector<pair<int, int>>();
                    vector<set<int>> *all_neighbors = new vector<set<int>>(n);

                    for (int i = 0; i < totalcount; i += 2)
                    {
                        all_edges->push_back({recvbuf[i], recvbuf[i + 1]});
                        (*all_neighbors)[recvbuf[i]].insert(recvbuf[i + 1]);
                        (*all_neighbors)[recvbuf[i + 1]].insert(recvbuf[i]);
                    }
                    // printf("hello2\n");

                    map<int, int> *visited = new map<int, int>();
                    queue<int> *trav = new queue<int>();
                    auto it = (*all_edges).begin();

                    while (true)
                    {
                        while (it != (*all_edges).end() && (*visited)[(*it).first] == 1 && (*visited)[(*it).second] == 1)
                            it++;
                        if (it == (*all_edges).end())
                            break;

                        int tmp_head = -1;
                        set<int> *grp_verts = new set<int>();
                        if ((*visited)[(*it).first] == 0)
                        {
                            (*trav).push((*it).first);
                            tmp_head = (*it).first;
                            (*visited)[(*it).first] = 1;
                        }
                        if ((*visited)[(*it).second] == 0)
                        {
                            (*trav).push((*it).second);
                            if (tmp_head == -1)
                                tmp_head = (*it).second;
                            (*visited)[(*it).second] = 1;
                        }
                        if (tmp_head == -1)
                            continue;
                        while ((*trav).size() > 0)
                        {
                            int i = (*trav).front();
                            (*trav).pop();
                            int node_rank = get_node_rank(i);

                            if (taskid == 2)
                            {
                                (*head)[i] = tmp_head;

                                // if (node_rank == 0)
                                //     (*head)[i] = tmp_head;
                                // else
                                // {
                                //     int node_head[2] = {i, tmp_head};
                                //     MPI_Send(&node_head, 2, MPI_INT, node_rank, 10, MPI_COMM_WORLD);
                                // }
                            }
                            grp_verts->insert(i);
                            (*visited)[i] = 1;
                            for (int j : (*all_neighbors)[i])
                            {
                                if ((*visited)[j] == 0)
                                {
                                    (*trav).push(j);
                                    (*visited)[j] = 1;
                                }
                            }
                        }
                        connected_comps->insert({tmp_head, *grp_verts});
                    }

                    if (taskid == 1)
                    {
                        if ((*connected_comps).size() == 0)
                            fout << 0 << endl;
                        else
                        {
                            fout << 1 << endl;
                            fout << (int)(*connected_comps).size() << endl;
                            for (auto grp : (*connected_comps))
                            {
                                int i = 0;
                                for (int v : grp.second)
                                {
                                    fout << v;
                                    if (i++ < (int)grp.second.size() - 1)
                                        fout << " ";
                                }
                                fout << endl;
                            }
                        }
                    }
                    else
                    {
                        // int node_head[2] = {-1, -1};
                        // for (int node_rank = 1; node_rank < size; node_rank++)
                        //     MPI_Send(&node_head, 2, MPI_INT, node_rank, 10, MPI_COMM_WORLD);

                        // /************** CALCULATING INFLUENCER SEQUENTIALLY **************

                        vector<set<int>> *conn_heads = new vector<set<int>>(n);
                        vector<int> *influencer = new vector<int>();
                        for (int vt = 0; vt < n; vt++)
                        {
                            if ((*head)[vt] != -1)
                                (*conn_heads)[vt].insert((*head)[vt]);
                            for (int i : (*neighbors)[vt])
                            {
                                if ((*head)[i] != -1)
                                    (*conn_heads)[vt].insert((*head)[i]);
                            }
                            if ((*conn_heads)[vt].size() >= task_p)
                                (*influencer).push_back(vt);
                        }
                        if (influencer->size() == 0)
                            fout << 0 << endl;
                        else
                        {
                            fout << influencer->size() << endl;
                            for (int in_vt : (*influencer))
                            {
                                fout << in_vt << " ";
                                if (verbose == 1)
                                {
                                    fout << "\n";
                                    for (int i : (*conn_heads)[in_vt])
                                        for (int j : (*connected_comps)[i])
                                            fout << j << " ";
                                    fout << endl;
                                }
                            }
                        }
                        // ************** CALCULATING INFLUENCER SEQUENTIALLY **************/
                    }
                }
                // printf("hello9\n");
                MPI_Barrier(MPI_COMM_WORLD);

                if (taskid == 2)
                {
                    //                     vector<set<int>> *conn_heads = new vector<set<int>>(num_nodes);
                    //                     if (rank == 0)
                    //                         conn_heads = new vector<set<int>>(n);

                    // #pragma omp parallel num_threads(2)
                    //                     {
                    //                         int tid = omp_get_thread_num();
                    //                         if (tid == 0)
                    //                         {
                    //                             for (int vert = first_node; vert < last_node; vert++)
                    //                             {
                    //                                 if ((*head)[vert] != -1)
                    //                                     (*conn_heads)[vert - first_node].insert((*head)[vert]);
                    //                                 for (int i : (*neighbors)[vert])
                    //                                 {
                    //                                     int node_rank = get_node_rank(i);
                    //                                     if (node_rank == rank && (*head)[i] != -1)
                    //                                         (*conn_heads)[vert - first_node].insert((*head)[i]);
                    //                                     else
                    //                                     {
                    //                                         int node_head[2] = {i, -1};
                    //                                         MPI_Send(&node_head, 1, MPI_INT, node_rank, 1, MPI_COMM_WORLD);
                    //                                         MPI_Recv(&node_head, 2, MPI_INT, node_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    //                                         assert(node_head[0] == i);
                    //                                         (*conn_heads)[vert - first_node].insert(node_head[1]);
                    //                                     }
                    //                                 }
                    //                             }
                    //                             int stop = -1;
                    //                             for (int oth_rank = 0; oth_rank < size; oth_rank++)
                    //                             {
                    //                                 if (oth_rank == rank)
                    //                                     continue;
                    //                                 MPI_Send(&stop, 1, MPI_INT, oth_rank, 1, MPI_COMM_WORLD);
                    //                             }
                    //                         }
                    //                         else
                    //                         {
                    //                             int node_head[2] = {-1, -1};
                    //                             MPI_Status status;
                    //                             int num_stop = size - 1;
                    //                             while (num_stop)
                    //                             {
                    //                                 MPI_Recv(&node_head, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
                    //                                 if (node_head[0] == -1)
                    //                                 {
                    //                                     num_stop--;
                    //                                     continue;
                    //                                 }
                    //                                 node_head[1] = (*head)[node_head[0] - first_node];
                    //                                 MPI_Send(&node_head, 2, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
                    //                             }
                    //                         }
                    //                     }
                    //                     if (rank != 0)
                    //                     {
                    //                         for (int vt = first_node; vt < last_node; vt++)
                    //                         {
                    //                             int vert = vt - first_node;
                    //                             if ((*conn_heads)[vert].size() < task_p)
                    //                                 continue;
                    //                             int sz = (*conn_heads)[vert].size();
                    //                             int conn_node[sz + 1];
                    //                             conn_node[0] = vt;
                    //                             int itr = 1;
                    //                             for (int ngbr : (*conn_heads)[vert])
                    //                                 conn_node[itr++] = ngbr;
                    //                             MPI_Send(&conn_node, 1 + sz, MPI_INT, 0, 9, MPI_COMM_WORLD);
                    //                         }
                    //                         int stop = -1;
                    //                         MPI_Send(&stop, 1, MPI_INT, 0, 9, MPI_COMM_WORLD);
                    //                     }
                    //                     else
                    //                     {
                    //                         int num_stop = size - 1;
                    //                         while (num_stop)
                    //                         {
                    //                             int tmp_conn[n];
                    //                             MPI_Recv(&tmp_conn, n, MPI_INT, MPI_ANY_SOURCE, 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    //                             if (tmp_conn[0] == -1)
                    //                             {
                    //                                 num_stop--;
                    //                                 continue;
                    //                             }
                    //                         }
                    //                     }
                }

                // end = chrono::system_clock::now();
                // elapsed_ms = end - start;
                // std::printf("Time to find %d traversal: %f milliseconds\n", k, 1000 * elapsed_ms.count());
            }
            else
            {
                vector<int> *head = new vector<int>(n, -1);

                // map<int, int> *head = new map<int, int>();
                map<int, set<int>> *headohead = new map<int, set<int>>();
                vector<vector<int>> *extr_vert = new vector<vector<int>>(size);

                // set<set<int>> *connected_comps = new set<set<int>>();

                map<int, set<int>> *connected_comps = new map<int, set<int>>();
                map<int, int> *visited = new map<int, int>();
                queue<int> *trav = new queue<int>();
                auto it = (*edges).begin();
                set<int> *all_heads;
                if (rank == 0)
                    all_heads = new set<int>();

                while (true)
                {
                    while (it != (*edges).end() && (*visited)[(*it).first] == 1 && (*visited)[(*it).second] == 1)
                        it++;
                    if (it == (*edges).end())
                        break;

                    int tmp_head = -1;
                    set<int> *grp_verts = new set<int>();
                    if ((*visited)[(*it).first] == 0)
                    {
                        (*trav).push((*it).first);
                        tmp_head = (*it).first;
                        (*visited)[(*it).first] = 1;
                    }
                    if ((*visited)[(*it).second] == 0)
                    {
                        (*trav).push((*it).second);
                        if (tmp_head == -1)
                            tmp_head = (*it).second;
                        (*visited)[(*it).second] = 1;
                    }
                    if (tmp_head == -1)
                        continue;

                    if (rank == 0)
                        all_heads->insert(tmp_head);

                    // printf("rank %d: tmp head %d\n", rank, tmp_head);
                    while ((*trav).size() > 0)
                    {
                        int i = (*trav).front();
                        (*trav).pop();
                        grp_verts->insert(i);
                        (*visited)[i] = 1;

                        int node_rank = get_node_rank(i);
                        (*head)[i] = tmp_head;
                        if (node_rank != rank)
                        {
                            (*extr_vert)[node_rank].push_back(i);
                            (*extr_vert)[node_rank].push_back(tmp_head);
                        }

                        for (int j : (*filtered_neighbors)[i])
                        {
                            if ((*visited)[j] == 0)
                            {
                                (*trav).push(j);
                                (*visited)[j] = 1;
                            }
                        }
                    }
                    connected_comps->insert({tmp_head, *grp_verts});
                    // connected_comps->insert(*grp_verts);
                }
                free(trav);
                free(visited);

                // printf("rank %d, k %d: cc size = %d\n", rank, k, connected_comps->size());
                // for (auto x : *connected_comps)
                //     printf("rank %d, k %d: head %d, size %d\n", rank, k, x.first, x.second.size());

                int num_stop = rank;
                vector<int> *tmp_conn = new vector<int>(2 * n);
                while (num_stop--)
                {
                    MPI_Status stat;
                    MPI_Recv(tmp_conn->data(), 2 * n, MPI_INT, MPI_ANY_SOURCE, 5, MPI_COMM_WORLD, &stat);
                    // printf("rank %d: received from %d\n", rank, stat.MPI_SOURCE);
                    if ((*tmp_conn)[0] != -1)
                    {
                        int cnt = 0;
                        MPI_Get_count(&stat, MPI_INT, &cnt);
                        for (int i = 0; i < cnt; i += 2)
                        {
                            if ((*head)[(*tmp_conn)[i]] == -1)
                                continue;
                            int hd = (*head)[(*tmp_conn)[i]];
                            // printf("head of %d is %d should be %d\n", (*tmp_conn)[i], hd, (*tmp_conn)[i + 1]);
                            // if (hd == (*tmp_conn)[i + 1])
                            //     continue;
                            // (*head)[hd] = (*tmp_conn)[i + 1];
                            (*headohead)[hd].insert((*tmp_conn)[i + 1]);
                        }
                    }
                }
                free(tmp_conn);

                for (int rk = rank + 1; rk < size; rk++)
                {
                    if ((*extr_vert)[rk].size() == 0)
                    {
                        int stop = -1;
                        MPI_Send(&stop, 1, MPI_INT, rk, 5, MPI_COMM_WORLD);
                    }
                    else
                    {
                        int sz = (*extr_vert)[rk].size();
                        vector<int> datum = (*extr_vert)[rk];
                        MPI_Send(datum.data(), sz, MPI_INT, rk, 5, MPI_COMM_WORLD);
                    }
                }

                // MPI_Barrier(MPI_COMM_WORLD);
                // printf("hello2 %d\n", rank);

                // printf("rank %d has size = %d for k = %d\n", rank, (int) connected_comps->size(),k);
                //  for(auto x:(*connected_comps))
                //  {
                //      for(int v:x)
                //      {
                //          printf("rank %d has vert = %d ",  rank,v);
                //      }
                //      printf("\n");
                //  }

                if (rank == 0)
                {
                    // map<int, int> mp_rep;
                    // map<int, set<int>> mp_comp;

                    // printf("rank 0 is here\n");

                    int num_stop = size - 1;
                    // for (int rk = 1; rk < size; rk++)
                    while (num_stop--)
                    {
                        while (true)
                        {
                            MPI_Status stat;
                            set<int> *recv_conn_comp = new set<int>();
                            int len;
                            // MPI_Recv(&len, 1, MPI_INT, rk, 6, MPI_COMM_WORLD, &stat);
                            MPI_Recv(&len, 1, MPI_INT, MPI_ANY_SOURCE, 6, MPI_COMM_WORLD, &stat);
                            if (len == -1)
                                break;
                            // if (len == -1) {
                            //     num_stop--;
                            //     continue;
                            // }
                            int vec[len];
                            int src = stat.MPI_SOURCE;
                            // MPI_Recv(vec, len, MPI_INT, rk, 7, MPI_COMM_WORLD, &stat);
                            MPI_Recv(vec, len, MPI_INT, src, 7, MPI_COMM_WORLD, &stat);

                            // int flag = 0;
                            // printf("Receive head of %d is %d\n", vec[1], vec[0]);
                            // printf("rank 0 got head %d from rank %d : %d\n", hd, stat.MPI_SOURCE, (*head)[hd]);

                            for (int id = 1; id < vec[0]; id++)
                            {
                                all_heads->insert(vec[id]);
                                (*headohead)[vec[vec[0]]].insert(vec[id]);
                                (*headohead)[vec[id]].insert(vec[vec[0]]);
                            }
                            all_heads->insert(vec[1]);
                            for (int ind = vec[0]; ind < len; ind++)
                            {
                                (*connected_comps)[vec[vec[0]]].insert(vec[ind]);
                                (*head)[vec[ind]] = vec[vec[0]];
                                // printf("hd = %d, vert = %d", vec[vec[0]], vec[ind]);
                            }

                            // if (head->find(vec[0]) == head->end())
                            //     (*head)[vec[0]] = vec[0];
                            // (*head)[vec[1]] = (*head)[vec[0]];
                            // int hd = (*head)[vec[0]];
                            // for (int ind = 1; ind < len; ind++)
                            //     (*connected_comps)[hd].insert(vec[ind]);
                            // for (auto x : (*connected_comps))
                            // printf("head = %d\n", x.first);
                        }
                    }

                    // printf("Finished stage!\n");

                    set<set<int>> *head_comps = new set<set<int>>();
                    map<int, int> *visited = new map<int, int>();
                    queue<int> *trav = new queue<int>();
                    auto it = all_heads->begin();

                    while (true)
                    {
                        while (it != all_heads->end() && (*visited)[*it] == 1)
                            it++;
                        if (it == all_heads->end())
                            break;

                        set<int> *grp_heads = new set<int>();
                        trav->push(*it);
                        (*visited)[*it] = 1;
                        // printf("rank %d: tmp head %d\n", rank, tmp_head);
                        while ((*trav).size() > 0)
                        {
                            int i = (*trav).front();
                            (*trav).pop();
                            grp_heads->insert(i);
                            (*visited)[i] = 1;

                            for (int j : (*headohead)[i])
                            {
                                if ((*visited)[j] == 0)
                                {
                                    trav->push(j);
                                    (*visited)[j] = 1;
                                }
                            }
                        }
                        head_comps->insert(*grp_heads);
                    }
                    free(visited);
                    free(trav);

                    if (taskid == 1)
                    {
                        map<int, int> *vis = new map<int, int>();
                        // printf("%d\n", connected_comps->size());
                        if ((*head_comps).size() == 0)
                            fout << 0 << endl;
                        else
                        {
                            fout << 1 << endl;
                            fout << (int)(*head_comps).size() << endl;
                            for (auto grp : (*head_comps))
                            {
                                // printf("k : %d, head grp: ", k);
                                for (int hd : grp)
                                {
                                    // printf("%d ", hd);
                                    // int sz = (*connected_comps)[hd].size() - 1;
                                    for (int v : (*connected_comps)[hd])
                                    {
                                        if ((*vis)[v] != 1)
                                        {
                                            fout << v << " ";
                                            (*vis)[v] = 1;
                                        }
                                        // sz--;
                                    }
                                }
                                // printf("\n");
                                fout << "\n";
                            }
                        }
                        free(vis);
                    }
                    else
                    {
                        // /************** CALCULATING INFLUENCER SEQUENTIALLY **************
                        map<int, int> *head_id = new map<int, int>();
                        int id = 0;
                        for (auto grp : (*head_comps))
                        {
                            for (int hd : grp)
                                (*head_id)[hd] = id;
                            id++;
                        }

                        vector<set<int>> *conn_head_ids = new vector<set<int>>(n);
                        vector<int> *influencer = new vector<int>();
                        for (int vt = 0; vt < n; vt++)
                        {
                            if ((*head)[vt] != -1)
                                (*conn_head_ids)[vt].insert((*head_id)[(*head)[vt]]);
                            for (int i : (*neighbors)[vt])
                                if ((*head)[i] != -1)
                                    (*conn_head_ids)[vt].insert((*head_id)[(*head)[i]]);
                            if ((*conn_head_ids)[vt].size() >= task_p)
                                (*influencer).push_back(vt);
                        }
                        if (influencer->size() == 0)
                            fout << 0 << endl;
                        else
                        {
                            fout << influencer->size() << endl;
                            int z = head_comps->size();
                            for (int in_vt : (*influencer))
                            {
                                if (verbose == 0)
                                    fout << in_vt << " ";
                                else
                                {
                                    map<int, int> *vis = new map<int, int>();
                                    fout << in_vt << "\n";
                                    auto itr = head_comps->begin();
                                    for (int i = 0; i < z; i++, itr++)
                                        if ((*conn_head_ids)[in_vt].find(i) != (*conn_head_ids)[in_vt].end())
                                            for (int h : (*itr))
                                                for (int v : (*connected_comps)[h])
                                                    if ((*vis)[v] != 1)
                                                    {
                                                        fout << v << " ";
                                                        (*vis)[v] = 1;
                                                    }
                                    // for (int i : (*conn_head_ids)[in_vt])
                                    //     for (int j : (*connected_comps)[i])
                                    //         fout << j << " ";
                                    fout << "\n";
                                    free(vis);
                                }
                            }
                        }
                        // ************** CALCULATING INFLUENCER SEQUENTIALLY **************/
                    }
                }
                else
                {
                    for (auto comp : *connected_comps)
                    {
                        auto it = comp.second.begin();
                        it++;
                        int hdsz = (*headohead)[comp.first].size() + 1;
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
                            vector<int> *tmp = new vector<int>(tmpsz);
                            int index = 0;
                            (*tmp)[index++] = hdsz;
                            for (auto x : (*headohead)[comp.first])
                                (*tmp)[index++] = x;
                            (*tmp)[index++] = comp.first;
                            // printf("1 size of %d is %d\n", comp.first, index);
                            for (; it != comp.second.end() && index < tmpsz; it++)
                                (*tmp)[index++] = *it;
                            // cout << "\n";
                            // printf("1 size of %d is %d\n", comp.first, tmpsz);
                            MPI_Send(&index, 1, MPI_INT, 0, 6, MPI_COMM_WORLD);
                            // printf("2 size of %d is %d\n", comp.first, tmp->size());
                            // printf("rank = %d sent %d size\n", rank, index);
                            // MPI_Buffer_attach(tmp.data(), index * sizeof(int));
                            MPI_Send(tmp->data(), index, MPI_INT, 0, 7, MPI_COMM_WORLD);
                            // printf("3 size of %d is %d\n", comp.first, index);
                            // printf("rank = %d sent %d size\n", rank, (int)comp.size());
                        }
                    }

                    //     for (auto comp : *connected_comps)
                    //     {
                    //         int hdsz = (*headohead)[comp.first].size() + 1;
                    //         vector<int> tmp(comp.second.size() + hdsz);
                    //         int index = 0;
                    //         tmp[index++] = hdsz;
                    //         for (auto x : (*headohead)[comp.first])
                    //             tmp[index++] = x;
                    //         printf("1 size of %d is %d\n", comp.first, index);
                    //         printf("1 size of %d is %d\n", comp.first, index + (int)comp.second.size());
                    //         for (auto v : comp.second)
                    //             tmp[index++] = v;
                    //         MPI_Send(&index, 1, MPI_INT, 0, 6, MPI_COMM_WORLD);
                    //         printf("2 size of %d is %d\n", comp.first, index);
                    //         // printf("rank = %d sent %d size\n", rank, index);
                    //         // MPI_Buffer_attach(tmp.data(), index * sizeof(int));
                    //         MPI_Send(tmp.data(), index, MPI_INT, 0, 7, MPI_COMM_WORLD);
                    //         printf("3 size of %d is %d\n", comp.first, index);
                    //         // printf("rank = %d sent %d size\n", rank, (int)comp.size());
                    //     }
                    int c = -1;
                    MPI_Send(&c, 1, MPI_INT, 0, 6, MPI_COMM_WORLD);
                    // printf("rank = %d sent -1\n", rank);
                }

                MPI_Barrier(MPI_COMM_WORLD);
                // }

                // set<set<int>> *connected_comps = new set<set<int>>();
                // map<int, int> *visited = new map<int, int>();
                // queue<int> *trav = new queue<int>();
                // auto it = (*edges).begin();

                // while (true)
                // {
                //     while (it != (*edges).end() && (*visited)[(*it).first] == 1 && (*visited)[(*it).second] == 1)
                //         it++;
                //     if (it == (*edges).end())
                //         break;

                //     set<int> *grp_verts = new set<int>();
                //     if ((*visited)[(*it).first] == 0)
                //     {
                //         (*trav).push((*it).first);
                //         (*visited)[(*it).first] = 1;
                //     }
                //     if ((*visited)[(*it).second] == 0)
                //     {
                //         (*trav).push((*it).second);
                //         (*visited)[(*it).second] = 1;
                //     }
                //     while ((*trav).size() > 0)
                //     {
                //         int i = (*trav).front();
                //         (*trav).pop();
                //         grp_verts->insert(i);
                //         (*visited)[i] = 1;
                //         for (int j : (*filtered_neighbors)[i])
                //         {
                //             if ((*visited)[j] == 0)
                //             {
                //                 (*trav).push(j);
                //                 (*visited)[j] = 1;
                //             }
                //         }
                //     }
                //     connected_comps->insert(*grp_verts);
                // }

                // printf("rank %d, k %d: cc size = %d\n", rank, k, connected_comps->size());

                // MPI_Barrier(MPI_COMM_WORLD);
                // printf("hello2 %d\n", rank);

                // printf("rank %d has size = %d for k = %d\n", rank, (int) connected_comps->size(),k);
                //  for(auto x:(*connected_comps))
                //  {
                //      for(int v:x)
                //      {
                //          printf("rank %d has vert = %d ",  rank,v);
                //      }
                //      printf("\n");
                //  }

                // if (rank == 0)
                // {

                //     map<int, int> mp_rep;
                //     map<int, set<int>> mp_comp;
                //     int rep = 0;
                //     for (auto comp : *connected_comps)
                //     {
                //         for (auto v : comp)
                //         {
                //             mp_rep[v] = rep;
                //         }
                //         mp_comp[rep] = comp;
                //         rep++;
                //     }

                //     // printf("rank 0 is here\n");

                //     int num_fin = 0;

                //     while (num_fin != size - 1)
                //     {
                //         MPI_Status stat;

                //         set<int> *recv_conn_comp = new set<int>();

                //         int len;

                //         MPI_Recv(&len, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &stat);

                //         // printf("rank 0 got1  %d\n",len);

                //         if (len == -1)
                //         {
                //             num_fin++;
                //             // printf("num_fin = %d\n",num_fin);
                //             continue;
                //         }
                //         int vec[len];
                //         MPI_Recv(vec, len, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &stat);

                //         // printf("rank 0 got2 %d\n",len);

                //         int flag = 0;

                //         set<set<int>> comps_to_be_del;
                //         set<int> comp_fin;
                //         set<int> visited;
                //         for (int in = 0; in < len; ++in)
                //         {
                //             int v = vec[in];
                //             if (mp_rep.find(v) != mp_rep.end())
                //             {
                //                 flag++;
                //                 auto comp = mp_comp[mp_rep[v]];

                //                 if (visited.find(mp_rep[v]) == visited.end())
                //                 {
                //                     comps_to_be_del.insert(comp);
                //                     comp_fin.insert(comp.begin(), comp.end());
                //                     visited.insert(mp_rep[v]);
                //                 }
                //             }
                //         }
                //         // printf("%d\n", flag);
                //         if (flag != 0)
                //         {
                //             for (int in = 0; in < len; ++in)
                //             {
                //                 comp_fin.insert(vec[in]);
                //             }
                //             for (auto comp : comps_to_be_del)
                //             {
                //                 connected_comps->erase(comp);
                //             }
                //             connected_comps->insert(comp_fin);
                //             for (auto v : comp_fin)
                //             {
                //                 mp_rep[v] = rep;
                //             }
                //             mp_comp[rep] = comp_fin;
                //             rep++;
                //         }
                //         else
                //         {
                //             set<int> recv_conn_comp;
                //             for (int in = 0; in < len; ++in)
                //             {
                //                 int v = vec[in];
                //                 mp_rep[v] = rep;
                //                 recv_conn_comp.insert(v);
                //             }
                //             mp_comp[rep] = recv_conn_comp;
                //             rep++;
                //             connected_comps->insert(recv_conn_comp);
                //         }
                //         // printf("%d\n", connected_comps->size());
                //     }
                //     // printf("%d\n", connected_comps->size());
                //     if ((*connected_comps).size() == 0)
                //         fout << 0 << endl;
                //     else
                //     {
                //         fout << 1 << endl;
                //         if (verbose == 1)
                //         {
                //             fout << (int)(*connected_comps).size() << endl;
                //             for (auto grp : (*connected_comps))
                //             {
                //                 int i = 0;
                //                 for (int v : grp)
                //                 {
                //                     fout << v;
                //                     if (i++ < (int)grp.size() - 1)
                //                         fout << " ";
                //                 }
                //                 fout << endl;
                //             }
                //         }
                //     }
                // }
                // else
                // {
                //     for (auto comp : *connected_comps)
                //     {
                //         vector<int> tmp(comp.size());
                //         int index = 0;
                //         for (auto v : comp)
                //         {
                //             tmp[index++] = v;
                //         }
                //         MPI_Ssend(&index, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                //         // printf("rank = %d sent %d size\n", rank, index);
                //         MPI_Ssend(tmp.data(), (int)comp.size(), MPI_INT, 0, 1, MPI_COMM_WORLD);
                //         // printf("rank = %d sent %d size\n", rank, (int)comp.size());
                //     }
                //     int c = -1;
                //     MPI_Ssend(&c, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                //     // printf("rank = %d sent -1\n", rank);
                // }

                // MPI_Barrier(MPI_COMM_WORLD);
            }
        }
    }
    if (rank == 0)
    {
        // fout.close();

        end_final = chrono::system_clock::now();
        chrono::duration<double> elapsed_ms = end_final - start_final;
        std::cout << "Time to end program: " << 1000 * elapsed_ms.count() << " milliseconds\n";
    }

    MPI_Finalize();
}
