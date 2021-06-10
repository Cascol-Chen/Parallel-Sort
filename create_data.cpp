#include<bits/stdc++.h>
#define INF 0x3f3f3f3f
using namespace std;
const int x=128;
fstream f;
template<typename T>
void wr(T r)
{
    f.write((char *)&r,sizeof(T));
}

template<typename T,typename ...Ts>
void wr(T r,Ts... args)
{
    f.write((char *)&r,sizeof(T));
    wr(args...);
}

int main()
{
    const int n=x*1000000;
    // f.open("test.txt",ios::out|ios::binary);
    // wr(3,6,7,1,0,2,4,5,8,9);
    f.open(to_string(x)+"M123.txt",ios::out|ios::binary);
    int res=INF;
    for(int i=1;i<=n;++i){
        int r=rand();
        f.write((char*)&r,sizeof(int));
    }
    f.close();
    // cout<<res<<"\n";
}