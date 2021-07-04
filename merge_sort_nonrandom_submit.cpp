/*
This is parallel merge_sort, the algorithm is as follow 
    1. assign a block to each process
    2. each process sort block independently
    3. find the k-1 split points in each process, and gather them up
    4. find k-1 split points among (k-1)*k split points
    5. split the array in each process into k segments according to split points and send them to the right process, each process will then have k segments
    6. each process run k-merge algorithm to merge k segments
    7. gather the final answer

author: Cascol-SCUT
link: https://github.com/Cascol-SCUT/Parallel-Sort
*/

#if 1
#pragma GCC optimize(1)
#pragma GCC optimize(2)
#pragma GCC optimize(3)
#pragma GCC optimize("Ofast")
#pragma GCC optimize("inline")
#pragma GCC optimize("-fgcse")
#pragma GCC optimize("-fgcse-lm")
#pragma GCC optimize("-fipa-sra")
#pragma GCC optimize("-ftree-pre")
#pragma GCC optimize("-ftree-vrp")
#pragma GCC optimize("-fpeephole2")
#pragma GCC optimize("-ffast-math")
#pragma GCC optimize("-fsched-spec")
#pragma GCC optimize("unroll-loops")
#pragma GCC optimize("-falign-jumps")
#pragma GCC optimize("-falign-loops")
#pragma GCC optimize("-falign-labels")
#pragma GCC optimize("-fdevirtualize")
#pragma GCC optimize("-fcaller-saves")
#pragma GCC optimize("-fcrossjumping")
#pragma GCC optimize("-fthread-jumps")
#pragma GCC optimize("-funroll-loops")
#pragma GCC optimize("-freorder-blocks")
#pragma GCC optimize("-fschedule-insns")
#pragma GCC optimize("inline-functions")
#pragma GCC optimize("-ftree-tail-merge")
#pragma GCC optimize("-fschedule-insns2")
#pragma GCC optimize("-fstrict-aliasing")
#pragma GCC optimize("-falign-functions")
#pragma GCC optimize("-fcse-follow-jumps")
#pragma GCC optimize("-fsched-interblock")
#pragma GCC optimize("-fpartial-inlining")
#pragma GCC optimize("no-stack-protector")
#pragma GCC optimize("-freorder-functions")
#pragma GCC optimize("-findirect-inlining")
#pragma GCC optimize("-fhoist-adjacent-loads")
#pragma GCC optimize("-frerun-cse-after-loop")
#pragma GCC optimize("inline-small-functions")
#pragma GCC optimize("-finline-small-functions")
#pragma GCC optimize("-ftree-switch-conversion")
#pragma GCC optimize("-foptimize-sibling-calls")
#pragma GCC optimize("-fexpensive-optimizations")
#pragma GCC optimize("inline-functions-called-once")
#pragma GCC optimize("-fdelete-null-pointer-checks")
#endif

#include <bits/stdc++.h>
#include <mpi.h>
using namespace std;

void break_point() //For Debug
{
    cout << "here\n";
    MPI_Finalize();
    exit(0);
}
struct node
{
    int v, id;
    bool operator<(const node &x) const
    {
        return v > x.v;
    }
};
const int mask = 0xffff, range = 1 << 16;
// 写内存连续，读内存跳跃（优化不大）
inline void radix_sort(vector<int> &a, bool f)
{
    vector<int> b(a.size());
    memcpy(&a[0], &b[0], sizeof(int) * a.size()); // 这句话是关键 优化很大 与缓存有关
    int lb[range], ub[range];
    memset(lb, 0, sizeof(lb));
    memset(ub, 0, sizeof(ub));
    for (const auto &it : a)
        lb[it & mask]++, ub[it >> 16]++;
    for (int i = 1; i < range; ++i)
        lb[i] += lb[i - 1];
    double st, ed;
    st = MPI_Wtime();
    for (const auto &it : a)
    {
        const int tmp = it & mask;
        b[lb[tmp] - 1] = it;
        lb[tmp]--;
    }
    ed = MPI_Wtime();
    // if(f)cout<<ed-st<<"s\n";
    lb[0] = 0;
    for (int i = 1; i < range; ++i)
        lb[i] = lb[i - 1] + ub[i - 1];
    for (const auto &it : b)
    {
        const int tmp = it >> 16;
        a[lb[tmp]++] = it;
    }
}
int main(int argc, char *argv[])
{
    int num_procs, rank;
    MPI_File fh;
    MPI_Status status;
    vector<int> recv, p, split_r, split_cnt;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    string file_name(argv[1]);
    int output_index = stoi(argv[2]);
    // file_name="test.txt";
    // Read in File
    if (MPI_File_open(MPI_COMM_WORLD, file_name.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh) != MPI_SUCCESS)
    {
        cout << "there is an error opening data\n";
        break_point();
    }

    MPI_Offset ofs;
    MPI_File_get_size(fh, &ofs);
    const int n = ofs / sizeof(int);
    const int number_per_process = n / num_procs;
    const int sz = (rank == num_procs - 1 ? (n - number_per_process * (num_procs - 1)) : number_per_process);
    int * const a = new int[sz];
    MPI_File_read_at(fh, rank * number_per_process * sizeof(int), &a[0], sz, MPI_INT, &status);

    MPI_File_close(&fh);

    // allocate all the required memory initially
    split_r.reserve(num_procs);
    split_cnt.reserve(num_procs);
    recv.resize(num_procs * (num_procs - 1));
    int cnt_per_process[num_procs], displ_per_process[num_procs], vec_info[num_procs][num_procs], displ[num_procs][num_procs];
    memset(cnt_per_process, 0, sizeof(cnt_per_process));
    // sort the block assigned

    // radix_sort(a,rank==0);
    sort(a, a + sz);
    // find how to split the array so as to make blocks with similar size
    const int per_block_size = sz / num_procs;
    p.resize(num_procs - 1);
    for (int i = per_block_size, j = 0; j < num_procs - 1; i += per_block_size, ++j)
        p[j] = a[i];

    MPI_Allgather(&p[0], num_procs - 1, MPI_INT, &recv[0], num_procs - 1, MPI_INT, MPI_COMM_WORLD);
    sort(recv.begin(), recv.end());
    for (int i = 0, j = 0; j < num_procs - 1; i += num_procs, ++j)
        p[j] = recv[i];

    int current_pos = 0, value_cnt = 0; // [l,r)
    for (int i = 0, len = p.size(); i < len; i += value_cnt)
    {
        value_cnt = 1;
        while (p[i + value_cnt] == p[i])
            ++value_cnt;
        const int pre_r = split_r.size() ? split_r.back() : 0;
        const auto it = lower_bound(a + pre_r, a + sz, p[i]);
        if (*it != p[i])
        { // Value doesn't exist in the array
            current_pos = it - a;
            for (int j = 0; j < value_cnt; ++j)
            {
                split_cnt.push_back(current_pos - (split_r.size() > 0 ? split_r.back() : 0));
                split_r.push_back(current_pos);
            }
        }
        else
        {
            const auto it2 = upper_bound(a, a + sz, p[i]);
            int per_num;
            const int total_number = it2 - a - pre_r;
            if (pre_r + total_number / (value_cnt + 1) < it - a)
            {
                per_num = (it2 - it) / (value_cnt + 1);
                current_pos = it - a + 1 + per_num;
            }
            else
            {
                per_num = total_number / (value_cnt + 1);
                current_pos = pre_r + per_num + 1;
            }
            for (int i = 0; i < value_cnt; ++i)
            {
                split_cnt.push_back(current_pos - (split_r.size() > 0 ? split_r.back() : 0));
                split_r.push_back(current_pos);
                current_pos += per_num;
            }
        }
    }
    split_cnt.push_back(sz - (split_r.size() > 0 ? split_r.back() : 0));
    split_r.push_back(sz);

    for (int i = 0; i < num_procs; ++i)
        MPI_Allgather(&split_cnt[i], 1, MPI_INT, &vec_info[i][0], 1, MPI_INT, MPI_COMM_WORLD);
    int total_cnt = 0;
    for (const auto &it : vec_info[rank])
        total_cnt += it;

    // show how many items are assigned to this process
    // cout<<rank<<": "<<total_cnt<<"\n";
    for (int i = 0; i < num_procs; ++i)
    {
        displ[i][0] = 0;
        for (int j = 1; j < num_procs; ++j)
        {
            displ[i][j] = displ[i][j - 1] + vec_info[i][j - 1];
        }
    }
    int * const recv_seg = new int[total_cnt];
    MPI_Request req[num_procs];
    int l = 0;
    for (int i = 0; i < num_procs; ++i)
    {
        MPI_Igatherv(&a[l], split_cnt[i], MPI_INT, &recv_seg[0], &vec_info[i][0], &displ[i][0], MPI_INT, i, MPI_COMM_WORLD, &req[i]);
        l = split_r[i];
    }
    MPI_Waitall(num_procs, req, MPI_STATUSES_IGNORE);

    // memory management
    delete[] a;
    int * const send = new int[total_cnt];

    int pl[num_procs], pr[num_procs];
    for (int i = 0; i < num_procs; ++i)
    {
        pl[i] = displ[rank][i];
        pr[i] = (i != num_procs - 1 ? displ[rank][i + 1] - 1 : total_cnt - 1);
    }
    // Merge K segments
    priority_queue<node> q;
    for (int i = 0; i < num_procs; ++i)
    {
        if (pl[i] <= pr[i])
        {
            q.push({recv_seg[pl[i]++], i});
        }
    }
    int idx = 0;
    while (!q.empty())
    {
        const auto tmp = q.top();
        const int id = tmp.id;
        q.pop();
        send[idx++] = tmp.v;
        while (pl[id] <= pr[id] && (q.empty() || q.top().v >= recv_seg[pl[id]]))
            send[idx++] = recv_seg[pl[id]++];
        if (pl[id] <= pr[id])
        {
            q.push({recv_seg[pl[id]++], id});
        }
    }

    // Gather final answer
    for (int i = 0; i < num_procs; ++i)
    {
        for (int j = 0; j < num_procs; ++j)
            cnt_per_process[i] += vec_info[i][j];
    }

    displ_per_process[0] = 0;
    for (int i = 1; i < num_procs; ++i)
        displ_per_process[i] = displ_per_process[i - 1] + cnt_per_process[i - 1];

    delete[] recv_seg;
    int * const ans = (!rank ? new int[n] : nullptr);
    MPI_Gatherv(&send[0], total_cnt, MPI_INT, &ans[0], &cnt_per_process[0], &displ_per_process[0], MPI_INT, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        cout<<ans[output_index+1]<<"\n";
    }
    MPI_Finalize();
}
/*
g++ -o .\merge_sort_nonrandom.exe .\merge_sort_nonrandom.cpp -l msmpi -L "C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64" -I "C:\Program Files (x86)\Microsoft SDKs\MPI\Include"
*/