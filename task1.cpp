#include <bits/stdc++.h>
#include <mpi.h>
#include <omp.h>

using namespace std;

typedef unsigned int uint;
typedef long long ll;

typedef struct
{
    int x;
    int y;
} Pair;

class TSQueue
{
public:
    void push2(int const &item1, int const &item2)
    {
        m.lock();
        m_queue.push_back(item1);
        m_queue.push_back(item2);
        m.unlock();
    }

    bool ispresent2(int &item1, int &item2)
    {
        m.lock();
        if (m_queue.size() < 2)
        {
            m.unlock();
            item1 = -1;
            item2 = -1;
            return false;
        }
        else
        {
            item1 = m_queue.front();
            m_queue.pop_front();
            item2 = m_queue.front();
            m_queue.push_front(item1);
            m.unlock();
            return true;
        }
    }
    bool pop2(int &item1, int &item2)
    {
        m.lock();
        if (m_queue.size() < 2)
        {
            m.unlock();
            item1 = -1;
            item2 = -1;
            return false;
        }
        else
        {
            item1 = m_queue.front();
            m_queue.pop_front();
            item2 = m_queue.front();
            m_queue.pop_front();
            m.unlock();
            return true;
        }
    }

    int size()
    {
        m.lock();
        int sz = m_queue.size();
        m.unlock();
        return sz;
    }

    void show_queue()
    {
        m.lock();
        for (int i = 0; i < m_queue.size(); i++)
        {
            printf("%d ", m_queue[i]);
        }
        printf("\n");
        m.unlock();
        return;
    }

private:
    deque<int> m_queue;
    mutex m;
};

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
    }
    if (inputpath.length() == 0 || headerpath.length() == 0)
    {
        std::cout << "Optional argument not provided\n";
        abort();
    }

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

    MPI_Datatype PairType;
    MPI_Type_contiguous(2, MPI_INT, &PairType);
    MPI_Type_commit(&PairType);

    /***************** GET INPUT FROM FILE *****************/

    chrono::time_point<std::chrono::system_clock> start, end, start_final, end_final;
    start = std::chrono::system_clock::now();
    start_final = std::chrono::system_clock::now();

    FILE *ptr = fopen(inputpath.c_str(), "rb"); // r for read, b for binary
    FILE *header = fopen(headerpath.c_str(), "rb");

    unsigned char buffer[4];

    fread(buffer, sizeof(buffer), 1, ptr); // read 8 bytes to our buffer
    int n = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);
    fread(buffer, sizeof(buffer), 1, ptr); // read 8 bytes to our buffer
    int m = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);

    set<pair<int, int>> *edges = new set<pair<int, int>>();

    int degi, node, p, offset, flag, node_bufj, nodej;

    if (size == 1)
    {
        for (int i = 0; i < n; ++i)
        {
            node = i, p = 0;
            fseek(header, 4 * node, SEEK_SET);
            fread(buffer, sizeof(buffer), 1, header);

            offset = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);
            fseek(ptr, offset + 4, SEEK_SET);
            fread(buffer, sizeof(buffer), 1, ptr);

            degi = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);

            unsigned char *buf = (unsigned char *)malloc(degi * sizeof(int));
            fread(buf, degi * sizeof(int), 1, ptr);

            for (int j = 0; j < degi; ++j)
            {
                nodej = (buf[p] & 0xFF) | ((buf[p + 1] & 0xFF) << 8) | ((buf[p + 2] & 0xFF) << 16) | ((buf[p + 3] & 0xFF) << 24);
                p += 4;

                if (node < nodej)
                {
                    if (edges->find({node, nodej}) == edges->end())
                    {
                        edges->insert({node, nodej});
                    }
                }
            }
            std::free(buf);
        }
    }
    else
    {
        if (rank == 0)
        {
            for (int i = 0; i < n; ++i)
            {
                node = i, p = 0;
                fseek(header, 4 * node, SEEK_SET);
                fread(buffer, sizeof(buffer), 1, header);

                offset = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);
                fseek(ptr, offset + 4, SEEK_SET);
                fread(buffer, sizeof(buffer), 1, ptr);

                degi = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);

                unsigned char *buf = (unsigned char *)malloc(degi * sizeof(int));
                fread(buf, degi * sizeof(int), 1, ptr);

                int flag = 0, node_bufj;
                for (int j = 0; j < degi; ++j)
                {
                    nodej = (buf[p] & 0xFF) | ((buf[p + 1] & 0xFF) << 8) | ((buf[p + 2] & 0xFF) << 16) | ((buf[p + 3] & 0xFF) << 24);
                    p += 4;

                    if (node < nodej)
                    {
                        if (edges->find({node, nodej}) == edges->end())
                        {
                            edges->insert({node, nodej});

                            if ((int)edges->size() == m / size)
                            {
                                flag = 1;
                                node_bufj = p;
                                break;
                            }
                        }
                    }
                }
                std::free(buf);
                if (flag == 1)
                {
                    MPI_Send(&node, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
                    MPI_Send(&node_bufj, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
                    break;
                }
            }
        }
        else if (rank == size - 1)
        {
            MPI_Recv(&node, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&node_bufj, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            p = node_bufj;
            for (int i = node; i < n; ++i)
            {
                node = i;
                fseek(header, 4 * node, SEEK_SET);
                fread(buffer, sizeof(buffer), 1, header);

                offset = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);
                fseek(ptr, offset + 4, SEEK_SET);
                fread(buffer, sizeof(buffer), 1, ptr);

                degi = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);

                unsigned char *buf = (unsigned char *)malloc(degi * sizeof(int));
                fread(buf, degi * sizeof(int), 1, ptr);

                flag = 0;

                for (int j = p / 4; j < degi; ++j)
                {
                    nodej = (buf[p] & 0xFF) | ((buf[p + 1] & 0xFF) << 8) | ((buf[p + 2] & 0xFF) << 16) | ((buf[p + 3] & 0xFF) << 24);
                    p += 4;

                    if (node < nodej)
                    {
                        if (edges->find({node, nodej}) == edges->end())
                        {
                            edges->insert({node, nodej});
                        }
                    }
                }
                std::free(buf);
                p = 0;
            }
        }
        else
        {
            MPI_Recv(&node, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&node_bufj, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            p = node_bufj;
            for (int i = node; i < n; ++i)
            {
                node = i;
                fseek(header, 4 * node, SEEK_SET);
                fread(buffer, sizeof(buffer), 1, header);

                offset = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);
                fseek(ptr, offset + 4, SEEK_SET);
                fread(buffer, sizeof(buffer), 1, ptr);

                degi = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);

                unsigned char *buf = (unsigned char *)malloc(degi * sizeof(int));
                fread(buf, degi * sizeof(int), 1, ptr);

                flag = 0;
                for (int j = p / 4; j < degi; ++j)
                {
                    nodej = (buf[p] & 0xFF) | ((buf[p + 1] & 0xFF) << 8) | ((buf[p + 2] & 0xFF) << 16) | ((buf[p + 3] & 0xFF) << 24);
                    p += 4;

                    if (node < nodej)
                    {
                        if (edges->find({node, nodej}) == edges->end())
                        {
                            edges->insert({node, nodej});
                            if ((int)edges->size() == m / size)
                            {
                                flag = 1;
                                node_bufj = p;
                                break;
                            }
                        }
                    }
                }
                std::free(buf);
                p = 0;
                if (flag == 1)
                {
                    MPI_Send(&node, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
                    MPI_Send(&node_bufj, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
                    break;
                }
            }
        }
    }
    vector<set<int>> *neighbors = new vector<set<int>>(n);
    // vector<bool> *visit = new vector<bool>(n, false);

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
                fread(buffer, sizeof(buffer), 1, header);

                offset = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);
                fseek(ptr, offset + 4, SEEK_SET);
                fread(buffer, 4, 1, ptr); // read 8 bytes to our buffer

                degi = (buffer[0] & 0xFF) | ((buffer[1] & 0xFF) << 8) | ((buffer[2] & 0xFF) << 16) | ((buffer[3] & 0xFF) << 24);

                unsigned char *buf = (unsigned char *)malloc(degi * sizeof(int));
                fread(buf, degi * sizeof(int), 1, ptr); // read degi * 4 bytes to our buffer

                p = 0;
                while (degi--)
                {
                    nodej = (buf[p] & 0xFF) | ((buf[p + 1] & 0xFF) << 8) | ((buf[p + 2] & 0xFF) << 16) | ((buf[p + 3] & 0xFF) << 24);
                    p += 4;
                    (*neighbors)[nd].insert(nodej);
                }
                std::free(buf);
            }
        }
    }
    std::fclose(ptr);
    std::fclose(header);

    MPI_Barrier(MPI_COMM_WORLD);

    end = chrono::system_clock::now();
    chrono::duration<double> elapsed_ms = end - start;
    std::cout << "Time to read input: " << 1000 * elapsed_ms.count() << " milliseconds\n";

    /***************** PRE-PROCESSING GRAPH *****************/

    start = chrono::system_clock::now();

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

    ofstream fout(outputpath, ios::out);

    // bool finished = false;

    for (int k = startk; k <= endk; k++)
    {
        TSQueue queue_send;
        TSQueue queue_recv;
        atomic<int> tag = 0;
        atomic<int> flag = 0;
        atomic<int> final_tag = 0;

        atomic<int> sent_minus1 = 0;

        auto it = (*edges).begin();

        while (it != edges->end())
        {
            pair<int, int> e = *it;
            int i = e.first;
            int j = e.second;
            if ((int)(triangles->at({i, j})).size() < k)
            {
                it = edges->erase(it);
                queue_send.push2(i, j);
                printf("rank %d -- taginc\n", rank);
                tag++;
            }
            else
                it++;
        }

#pragma omp parallel num_threads(3)
        {
            int tid = omp_get_thread_num();

            if (tid == 0) // send thread
            {
                while (flag == 0)
                {
                    int e1, e2;

                    if (queue_send.pop2(e1, e2))
                    {
                        Pair pii;
                        pii.x = e1;
                        pii.y = e2;

                        for (int i = 0; i < size; ++i)
                        {
                            {
                                queue_send.show_queue();
                                printf("rank %d sent {%d,%d} to %d with tag %d\n", rank, e1, e2, i, tag.load());
                                final_tag = tag.load();
                                MPI_Send(&pii, 1, PairType, i, tag, MPI_COMM_WORLD);
                            }
                        }
                        sent_minus1 = 0;
                    }
                    else
                    {
                        if (sent_minus1 == 0 && queue_recv.size() == 0)
                        {
                            Pair pii;
                            pii.x = -1;
                            pii.y = -1;
                            printf("tid 0 sent -1\n");

                            for (int i = 0; i < size; ++i)
                            {
                                // if (tag > final_tag)
                                {
                                    printf("rank %d  sent -1 to %d with tag %d\n", rank, i, tag.load());
                                    final_tag = tag.load();
                                    MPI_Send(&pii, 1, PairType, i, tag, MPI_COMM_WORLD);
                                }
                            }
                            sent_minus1 = 1;
                        }
                    }
                }
            }

            if (tid == 1)
            {
                MPI_Status status;
                vector<int> latest_tags(size, -1);

                while (flag == 0)
                {
                    Pair pii;
                    printf("abt to rcv %d\n", rank);
                    MPI_Recv(&pii, 1, PairType, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                    printf("rank %d recvd {%d,%d} with tag %d from src %d\n", rank, pii.x, pii.y, status.MPI_TAG, status.MPI_SOURCE);

                    if (pii.x == -1)
                    {
                        assert(status.MPI_TAG >= latest_tags[status.MPI_SOURCE]);
                        latest_tags[status.MPI_SOURCE] = status.MPI_TAG;

                        bool stop = true;
                        for (int i = 0; i < size; ++i)
                        {
                            printf("tag in me = %d\n", tag.load());
                            if (latest_tags[i] != tag)
                            {
                                stop = false;
                                break;
                            }
                        }
                        printf("helo\n");

                        if (stop && queue_send.size() == 0 && queue_recv.size() == 0)
                        {
                            flag = 1;
                            printf("tag = %d\n", tag.load());
                            printf("rank %d broken\n", rank);
                            break;
                        }
                    }
                    else
                    {
                        queue_recv.push2(pii.x, pii.y);
                        if (status.MPI_SOURCE != rank)
                            tag++;
                    }
                }
            }

            if (tid == 2)
            {
                while (flag == 0)
                {
                    int i, j;
                    if (queue_recv.ispresent2(i, j))
                    {
                        for (int p : (*triangles)[{i, j}])
                        {
                            int w = min(i, p);
                            int x = max(i, p);
                            int y = min(j, p);
                            int z = max(j, p);
                            if ((*triangles)[{w, x}].find(j) != (*triangles)[{w, x}].end())
                                (*triangles)[{w, x}].erase(j);

                            if (edges->find({w, x}) != edges->end() && (int)(triangles->at({w, x})).size() < k)
                            {
                                edges->erase({w, x});
                                queue_send.push2(w, x);
                                tag++;
                            }

                            if ((*triangles)[{y, z}].find(i) != (*triangles)[{y, z}].end())
                                (*triangles)[{y, z}].erase(i);

                            if (edges->find({y, z}) != edges->end() && (int)(triangles->at({y, z})).size() < k)
                            {
                                edges->erase({y, z});
                                queue_send.push2(y, z);
                                tag++;
                            }
                        }

                        queue_recv.pop2(i, j);
                        if (tag.load() > final_tag)
                            sent_minus1 = 0;
                        (*neighbors)[i].erase(j);
                        (*neighbors)[j].erase(i);
                    }
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        end = chrono::system_clock::now();
        elapsed_ms = end - start;
        std::printf("Time to find %d-truss: %f milliseconds\n", k, 1000 * elapsed_ms.count());

        if (verbose == 0)
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
            start = chrono::system_clock::now();

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

            if (rank != 0)
            {
                MPI_Gather(&sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
                // printf("rank %d has %d count\n", rank, sendcount);
                MPI_Gatherv(edge_to_vert->data(), sendcount, MPI_INT, recvbuf, recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
                // printf("rank %d has %d count\n", rank, sendcount);
            }
            else
            {
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
                MPI_Gatherv(edge_to_vert->data(), sendcount, MPI_INT, recvbuf, recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

                vector<pair<int, int>> *all_edges = new vector<pair<int, int>>();
                vector<set<int>> *all_neighbors = new vector<set<int>>(n);

                for (int i = 0; i < totalcount; i += 2)
                {
                    all_edges->push_back({recvbuf[i], recvbuf[i + 1]});
                    (*all_neighbors)[recvbuf[i]].insert(recvbuf[i + 1]);
                    (*all_neighbors)[recvbuf[i + 1]].insert(recvbuf[i]);
                }
                // printf("hello2\n");

                set<set<int>> *connected_comps = new set<set<int>>();
                map<int, int> *visited = new map<int, int>();
                queue<int> *trav = new queue<int>();
                auto it = (*all_edges).begin();

                while (true)
                {
                    while (it != (*all_edges).end() && (*visited)[(*it).first] == 1 && (*visited)[(*it).second] == 1)
                        it++;
                    if (it == (*all_edges).end())
                        break;

                    set<int> *grp_verts = new set<int>();
                    if ((*visited)[(*it).first] == 0)
                    {
                        (*trav).push((*it).first);
                        (*visited)[(*it).first] = 1;
                    }
                    if ((*visited)[(*it).second] == 0)
                    {
                        (*trav).push((*it).second);
                        (*visited)[(*it).second] = 1;
                    }
                    while ((*trav).size() > 0)
                    {
                        int i = (*trav).front();
                        (*trav).pop();
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
                    connected_comps->insert(*grp_verts);
                }

                if ((*connected_comps).size() == 0)
                    fout << 0 << endl;
                else
                {
                    fout << 1 << endl;
                    if (verbose == 1)
                    {
                        fout << (int)(*connected_comps).size() << endl;
                        for (auto grp : (*connected_comps))
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
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    if (rank == 0)
    {
        // fout.close();

        end_final = chrono::system_clock::now();
        elapsed_ms = end_final - start_final;
        std::cout << "Time to end program: " << 1000 * elapsed_ms.count() << " milliseconds\n";
    }

    MPI_Finalize();
}
