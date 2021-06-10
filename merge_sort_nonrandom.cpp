/*
This is parallel merge_sort, the algorithm is as follow 
    1. assign a block to each process
    2. each process sort block independently
    3. find the k-1 split points in each process, and gather them up
    4. find k-1 split points among (k-1)*k split points
    5. split the array in each process into k segments according to split points and send them to the right process, each process will then have k segments
    6. each process run k-merge algorithm to merge k segments
    7. gather the final answer
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

#include<bits/stdc++.h>
#include<mpi.h>
#define max_thread 12
using namespace std;

const int x=256;
const int n=x*1e6;
// const int n=10;

void break_point() //调试使用函数
{
    cout<<"here\n";
    MPI_Finalize();
    exit(0);
}
struct node{
    int v,id;
    bool operator <(const node &x) const{
        return v>x.v;
    }
};
int main(int argc,char *argv[])
{
    int num_procs,rank;
    MPI_File fh;
    MPI_Status status;
    vector<int> a,vec_cnt,recv,p,split_r,split_cnt;
    vector<vector<int>> vec,vec_info,displ;
    double cst,ced,ed2,ed3,ed4,ed5,ed6,ed7,ed8,ed9,ed10;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&num_procs);

    const int number_per_process=n/num_procs;
    int sz=number_per_process;
    string file_name=to_string(x)+"M.txt";
    // file_name="test.txt";
    
    // Read in File
    if(MPI_File_open(MPI_COMM_WORLD,file_name.c_str(),MPI_MODE_RDONLY,MPI_INFO_NULL,&fh)!=MPI_SUCCESS){
        cout<<"there is an error opening data\n";
        break_point();
    }
    if(rank==num_procs-1) {
        sz=number_per_process+(n-number_per_process*num_procs);
        a.resize(sz);
        MPI_File_read_at_all(fh,rank*number_per_process*sizeof(int),&a[0],sz,MPI_INT,&status);
    }
    else{
        a.resize(number_per_process);
        MPI_File_read_at_all(fh,rank*number_per_process*sizeof(int),&a[0],number_per_process,MPI_INT,&status);
    }
    MPI_File_close(&fh);

    // allocate all the required memory initially
    split_r.reserve(num_procs);split_cnt.reserve(num_procs);
    recv.resize(num_procs*(num_procs-1));
    vec_info=vector<vector<int>>(num_procs,vector<int>(num_procs));
    displ=vector<vector<int>>(num_procs,vector<int>(num_procs));
    vector<int> cnt_per_process(num_procs),displ_per_process(num_procs),ans;
    if(rank==0) ans.resize(n);
    
    cst=MPI_Wtime();
    // sort the block assigned
    sort(a.begin(),a.end());
    ed2=MPI_Wtime();

    // find how to split the block so as to make them have similar size
    int per_block_size=sz/num_procs; p.resize(num_procs-1);
    for(int i=per_block_size,j=0;j<num_procs-1;i+=per_block_size,++j) p[j]=a[i];

    MPI_Allgather(&p[0],num_procs-1,MPI_INT,&recv[0],num_procs-1,MPI_INT,MPI_COMM_WORLD);
    sort(recv.begin(),recv.end());
    for(int i=0,j=0;j<num_procs-1;i+=num_procs,++j) p[j]=recv[i];
    ed4=MPI_Wtime();

    int current_pos=0, value_cnt=0; // [l,r)

    for(int i=0,len=p.size();i<len;i+=value_cnt)
    {
        value_cnt=1;
        while(p[i+value_cnt]==p[i]) ++value_cnt;
        auto it=lower_bound(a.begin(),a.end(),p[i]);
        if(*it!=p[i]) {  //没有找到
            current_pos=it-a.begin();
            for(int j=0;j<value_cnt;++j) {
                split_cnt.push_back(current_pos-(split_r.size()>0?split_r.back():0));
                split_r.push_back(current_pos);
            }
        }
        else{
            auto it2=upper_bound(a.begin(),a.end(),p[i]);
            int per_num;
            const int pre_r=split_r.size()?split_r.back():0;
            const int total_number=it2-a.begin()-pre_r;
            if(pre_r+total_number/(value_cnt+1)<it-a.begin()){
            // if(true){
                per_num=(it2-it)/(value_cnt+1);
                current_pos=it-a.begin()+1+per_num;
            }
            else{
                // assert(false);
                per_num=total_number/(value_cnt+1);
                current_pos=pre_r+per_num+1;
            }
            for(int i=0;i<value_cnt;++i){
                split_cnt.push_back(current_pos-(split_r.size()>0?split_r.back():0));
                split_r.push_back(current_pos);
                current_pos+=per_num;
            }
        }
    }
    split_cnt.push_back(sz-(split_r.size()>0?split_r.back():0));
    split_r.push_back(sz);
    ed5=MPI_Wtime();

    for(int i=0;i<num_procs;++i) MPI_Allgather(&split_cnt[i],1,MPI_INT,&vec_info[i][0],1,MPI_INT,MPI_COMM_WORLD);
    int total_cnt=0;
    for(const auto& it:vec_info[rank]) total_cnt+=it;
    cout<<rank<<": "<<total_cnt<<"\n";
    for(int i=1;i<num_procs;++i){
        for(int j=0;j<num_procs;++j){
            displ[j][i]=displ[j][i-1]+vec_info[j][i-1];
        }
    }
    vector<int> recv_seg(total_cnt);
    MPI_Request req[num_procs];
    int l=0;
    for(int i=0;i<num_procs;++i){
        MPI_Igatherv(&a[l],split_cnt[i],MPI_INT,&recv_seg[0],&vec_info[i][0],&displ[i][0],MPI_INT,i,MPI_COMM_WORLD,&req[i]);
        l=split_r[i];
    }
    MPI_Waitall(num_procs,req,MPI_STATUSES_IGNORE);

    // 先释放a的内存然后再申请内存
    vector<int>().swap(a); //释放a数组内存
    vector<int>().swap(recv);
    vector<int> send(total_cnt);
    ed6=MPI_Wtime();
    
    int pl[num_procs],pr[num_procs];
    for(int i=0;i<num_procs;++i){
        pl[i]=displ[rank][i];
        if(i+1<num_procs) pr[i]=displ[rank][i+1]-1;
        else pr[i]=total_cnt-1;
    }
    // merge
    priority_queue<node> q;
    for(int i=0;i<num_procs;++i){
        if(pl[i]<=pr[i]){
            q.push({recv_seg[pl[i]++],i});
        }
    }
    int idx=0;
    while(!q.empty())
    {
        auto tmp=q.top();q.pop();
        send[idx++]=tmp.v;
        while(pl[tmp.id]<=pr[tmp.id]&&(q.empty()||q.top().v>=recv_seg[pl[tmp.id]]) ) send[idx++]=recv_seg[pl[tmp.id]++]; 
        if(pl[tmp.id]<=pr[tmp.id]){
            q.push({recv_seg[pl[tmp.id]++],tmp.id});
        }
    }
    ed7=MPI_Wtime();
    // break_point();
    for(int i=0;i<num_procs;++i){
        for(int j=0;j<num_procs;++j) cnt_per_process[i]+=vec_info[i][j];
    }
    for(int i=1;i<num_procs;++i) displ_per_process[i]=displ_per_process[i-1]+cnt_per_process[i-1];
    MPI_Gatherv(&recv_seg[0],recv_seg.size(),MPI_INT,&ans[0],&cnt_per_process[0],&displ_per_process[0],MPI_INT,0,MPI_COMM_WORLD);
    ed8=MPI_Wtime();

    if(rank==0){
        ced=MPI_Wtime();
        // for(int i=0;i<1100;++i) cout<<ans[i]<<" "
        fstream cou;cou.open("record.txt",ios::app|ios::out);
        cout<<"Sort time: "<<ed2-cst<<"s\n";
        // cout<<"Get split value: "<<ed4-ed3<<"s\n";
        // cout<<"Get seg: "<<ed5-ed4<<"s\n";
        cout<<"Gather seg: "<<ed6-ed5<<"s\n";
        cout<<"Merge seg: "<<ed7-ed6<<"s\n";
        cout<<"Merge and Sort: "<<ed7-ed6+ed2-cst<<"s\n";
        cout<<"Final gather: "<<ed8-ed7<<"s\n";
        cout<<"Parallelize with -n "<<num_procs<<": "<<ced-cst<<"\n";

        cou<<"with -n "<<num_procs<<"\n";
        cou<<"Sort time: "<<ed2-cst<<"s\n";
        cou<<"Gather seg: "<<ed6-ed5<<"s\n";
        cou<<"Merge seg: "<<ed7-ed6<<"s\n";
        cou<<"Merge and Sort: "<<ed7-ed6+ed2-cst<<"s\n";
        cou<<"Final gather: "<<ed8-ed7<<"s\n";
        cou<<"Parallelize with -n "<<num_procs<<": "<<ced-cst<<"\n";
        cou<<"-----------------------\n";
    }
    MPI_Finalize();
}
/*
g++ -o .\merge_sort_nonrandom.exe .\merge_sort_nonrandom.cpp -l msmpi -L "C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64" -I "C:\Program Files (x86)\Microsoft SDKs\MPI\Include"
*/

/*
1
Sort time: 5.81909s
Merge seg: 2.81987s
Merge and Sort: 8.63896s

2
Sort time: 2.8936s
Merge seg: 2.00338s
Merge and Sort: 4.89699s

4
Sort time: 1.65492s
Merge seg: 1.4322s
Merge and Sort: 3.08712s

6
Sort time: 1.16097s
Merge seg: 1.07781s
Merge and Sort: 2.23877s

8
Sort time: 1.01886s
Merge seg: 1.1054s
Merge and Sort: 2.12425s

12
Sort time: 0.810835s
Merge seg: 1.07088s
Merge and Sort: 1.88171s
*/