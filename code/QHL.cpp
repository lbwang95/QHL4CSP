#include<bits/stdc++.h>
#include<unordered_map>
using namespace std;
typedef pair<double, double> DD;
typedef pair<double, int> DI;
typedef pair<int, vector<DD>> IV;
typedef pair<int, vector<int>> IVI;
typedef pair<int, int> II;
const int MAX_V = 435676;
int N, rootpa;
long long npc, nhoplink;
double optw = DBL_MAX;
int treeheight = 0, treewidth = 0;

typedef struct node{
    int level;
    int ran;
    int parent;
    vector<int> children;
    vector<int> X;
};
node T[MAX_V];
vector<IV> L[MAX_V];
unordered_map<int, int> TX[MAX_V];
vector<int> ancarray[MAX_V];
int root = -1;

void Sky(vector<DD> &P1,vector<DD> &P2, vector<DD> &res){
    //return the res contains the skyline paths of joining P1 and P2
    if(P1.size()==0)
        for (int i = 0; i < P2.size();i++)
            res.push_back(P2[i]);
    if(P2.size()==0)
        for (int i = 0; i < P1.size();i++)
            res.push_back(P1[i]);
    for (int i = 0; i < P1.size();i++)
        for (int j = 0; j < P2.size();j++){
            res.push_back(DD(P1[i].first + P2[j].first, P1[i].second + P2[j].second));
        }
    sort(res.begin(), res.end());
    double prev = DBL_MAX;
    vector<DD> tmp;
    for (int i = 0; i < res.size();i++){
        if(res[i].second<prev){
            tmp.push_back(res[i]);
            prev = res[i].second;
        }
    }
    res = tmp;
}

void load(string s){
    ifstream inf;
    string filename = string("../data/") + s + string("/index");
    inf.open(filename.c_str(), ios::binary);
    //fscanf(fp, "%d", &N);
    inf.read(reinterpret_cast<char *>(&N), sizeof(int));
    for (int k = 0; k < N; k++){
        int v, nx, lenl;
        //fscanf(fp, "%d", &v);
        //fscanf(fp, "%d%d%d", &T[v].parent, &nx, &lenl);
        inf.read(reinterpret_cast<char *>(&v), sizeof(int));
        inf.read(reinterpret_cast<char *>(&T[v].parent), sizeof(int));
        inf.read(reinterpret_cast<char *>(&nx), sizeof(int));
        inf.read(reinterpret_cast<char *>(&lenl), sizeof(int));
        if (k == 0){
            root = v;
            rootpa = T[v].parent;
        }
        for (int i = 0; i < nx; i++){
            int a;
            //fscanf(fp, "%d", &a);
            inf.read(reinterpret_cast<char *>(&a), sizeof(int));
            T[v].X.push_back(a);
            TX[v].insert(II(a, 1));
        }
        for (int i = 0; i < lenl; i++){
            int lend, a, b, d;
            //fscanf(fp, "%d%d", &a, &lend);
            inf.read(reinterpret_cast<char *>(&a), sizeof(int));
            inf.read(reinterpret_cast<char *>(&lend), sizeof(int));
            vector<DD> tmp;
            L[v].push_back(IV(a, tmp));
            if (TX[v].count(a) != 0)
                ancarray[v].push_back(i);
            for (int j = 0; j < lend;j++){
                //fscanf(fp, "%lf%lf", &b, &d);
                inf.read(reinterpret_cast<char *>(&b), sizeof(int));
                inf.read(reinterpret_cast<char *>(&d), sizeof(int));
                L[v][L[v].size() - 1].second.push_back(DD(b, d));
            }
        }
    }
    inf.close();
    //fclose(fp);
}

void concat(vector<DD> &P1,vector<DD> &P2, vector<DD> &res){
    if(P1.size()==0)
        for (int i = 0; i < P2.size();i++)
            res.push_back(P2[i]);
    if(P2.size()==0)
        for (int i = 0; i < P1.size();i++)
            res.push_back(P1[i]);
    for (int i = 0; i < P1.size();i++)
        for (int j = 0; j < P2.size();j++){
            res.push_back(DD(P1[i].first + P2[j].first, P1[i].second + P2[j].second));
        }
}

int curcs;
bool cmpp(const II &a,const II &b){
    return L[curcs][a.first].second[0].first < L[curcs][b.first].second[0].first;
}

void Cub_h(int v,int u,int h,vector<DD> &P_, vector<DD> &P__, double &Ch){
    //setting Cub_h for fixed u (Alg.5)
    int j = 0;
    for (int i = 0; i < P_.size();i++){
        while(j<P__.size()){
            if(P_[i].first==P__[j].first&&P_[i].second==P_[j].second)
                break;
            else if(P_[i].first<P__[j].first){
                Ch = max(Ch, P_[i].first);
                return;
            }
            else
                j++;
        }
        if(j==P__.size()){
            Ch = max(Ch, P_[i].first);
            return;
        }
    }
    Ch = DBL_MAX;
    return;
}
unordered_map<int, double> relationC[MAX_V];//the speedup technique
unordered_map<int,vector<DI>> Cub[MAX_V];//following ancarray positions -1
int totpcs = 0, totvisitpcs = 0;
long long qhlindexsize;
void pruningConditions(int cs, int v){
    //Alg. 6
    //printf("cs%d v%d\n", cs + 1, v + 1);
    if(Cub[v].count(cs)==1)
        return;
    totpcs += 1;
    vector<II> _anc;
    curcs = cs;
    for (int i = 0; i < ancarray[cs].size() - 1; i++){//minus cs
        Cub[v][cs].push_back(DI(0, ancarray[cs][i]));
        _anc.push_back(II(ancarray[cs][i], i));
    }
    sort(_anc.begin(), _anc.end(), cmpp);
    int lena = _anc.size();
    default_random_engine gen(time(NULL));
    for (int j = 1; j < lena; j++){
        uniform_int_distribution<int> iu(0, j-1);
        int i = iu(gen);
        vector<DD> res;
        int levu = _anc[i].first, levh = _anc[j].first;
        int u = L[cs][levu].first, h = L[cs][levh].first;
        if (relationC[v].count(levu * treeheight + levh) != 0)
        {
            Cub[v][cs][_anc[j].second].first = max(Cub[v][cs][_anc[j].second].first, relationC[v][levu * treeheight + levh]);
        }
        else{
            Sky(L[v][levu].second, (levu>levh)?L[u][levh].second:L[h][levu].second, res);
            Cub_h(v,u,h,L[v][levh].second, res, Cub[v][cs][_anc[j].second].first);
            relationC[v][levu * treeheight + levh] = Cub[v][cs][_anc[j].second].first;
            //printf("|%d %d %d %d %f|", levu, levh, u + 1, h + 1, Cub[v][cs][_anc[j].second]);
        }
    }
    qhlindexsize += Cub[v][cs].size();
}

void QHLindex(string prefix){//building all pruning conditions Sec. 4.2
    FILE *fpq = fopen((prefix + string("random")).c_str(), "r");
    int s = -1, t = 0;
    double C;
    int iter = 0;
    default_random_engine gen(time(NULL));
    uniform_int_distribution<int> st(0, N - 1);
    while (~fscanf(fpq, "%d%d%lf", &s, &t, &C)){
        iter += 1;
        if(iter%10000==0)
            printf("%d\n", iter);
        vector<int> ancs,anct;
        int u1 = s, l = -1;
        while(u1!=MAX_V){
            ancs.push_back(u1);
            u1 = T[u1].parent;
        }
        u1 = t;
        while(u1!=MAX_V){
            anct.push_back(u1);
            u1 = T[u1].parent;
        }
        int i = ancs.size() - 1, j = anct.size() - 1, k = -1;
        while (i != -1 && j != -1){
            if (ancs[i] == anct[j]){
                i--,j--,k++;
            }
            else
                break;
        }
        if(i==-1)
            l = ancs[0];
        else if(j==-1)
            l = anct[0];
        else
            l = ancs[i + 1];
        int ind;
        if (l == s||l == t){
            continue;
        }
        else{
            int cs = ancs[i], ct = anct[j];
            pruningConditions(cs, s);
            pruningConditions(cs, t);
            pruningConditions(ct, s);
            pruningConditions(ct, t);
        }
    }
    printf("Total number of pruning conditions %d\n", totpcs);
}
long long qhltc;
void pathConcatenationC(vector<DD> &Psh,vector<DD> &Pht, double &C){
    int i = 0, j = Pht.size() - 1;
    while (i != Psh.size() && j != -1){
        if (Psh[i].first + Pht[j].first <= C){
            qhltc++;
            if (Psh[i].second + Pht[j].second < optw){
                optw = Psh[i].second + Pht[j].second;
            }
            i++;
        }
        else
            j--;
    }
}
int hitpc;
void usePC(bool initflag, int c, int v, int s, int t, double C, vector<int> &H, int &retc, int &mtc){
    if(Cub[v].count(c) != 0){
        hitpc++;
        double tc = 0;
        vector<int> H1;
        for (int i = 0; i < Cub[v][c].size();i++){
            if(C>=Cub[v][c][i].first){
                int levh = ancarray[c][i];//forbid sort since ancarray follows Cub
                tc+=L[s][levh].second.size() + L[t][levh].second.size();
                H1.push_back(ancarray[c][i]);
            }
        }
        if(tc<mtc){
            H = H1;
            mtc = tc;
            retc = c;
        }
    }
    else if(initflag){
        mtc = 0;
        retc = c;
        for (int i = 0; i < ancarray[c].size() - 1; i++){//minus cs
            H.push_back(ancarray[c][i]);
            int levh = ancarray[c][i];
            mtc+=L[s][levh].second.size() + L[t][levh].second.size();
        }
    }
    totvisitpcs++;
}
int qhlhopsize;
void QHL(int s, int t, double C){
    optw = DBL_MAX;
    if (s == t)
        return;
    s--,t--;
    vector<int> ancs,anct;
    int u1 = s, l = -1;
    while(u1!=MAX_V){
        ancs.push_back(u1);
        u1 = T[u1].parent;
    }
    u1 = t;
    while(u1!=MAX_V){
        anct.push_back(u1);
        u1 = T[u1].parent;
    }
    int i = ancs.size() - 1, j = anct.size() - 1, k = -1;
    unordered_map<int, int> inds, indt;
    while (i != -1 && j != -1){
        if (ancs[i] == anct[j]){
            inds[ancs[i]] = i;
            indt[anct[j]] = j;
            i--, j--, k++;
        }
        else
            break;
    }
    if(i==-1)
        l = ancs[0];
    else if(j==-1)
        l = anct[0];
    else
        l = ancs[i + 1];
    int ind;
    if (l == s){
        IV iv = L[t][ancs.size() - 1];
        ind = upper_bound(iv.second.begin(), iv.second.end(), DD(C, DBL_MAX)) - iv.second.begin();
        if(ind!=0)
            optw = iv.second[ind - 1].second;
    }
    else if (l == t){
        IV iv = L[s][anct.size() - 1];
        ind=upper_bound(iv.second.begin(), iv.second.end(), DD(C, DBL_MAX))-iv.second.begin();
        if(ind!=0)
            optw = iv.second[ind - 1].second;
    }
    else{
        int cs = ancs[i], ct = anct[j];
        int c, mintc = INT_MAX;
        vector<int> H;
        usePC(1, cs, s, s, t, C, H, c, mintc);
        usePC(0, cs, t, s, t, C, H, c, mintc);
        usePC(0, ct, s, s, t, C, H, c, mintc);
        usePC(0, ct, t, s, t, C, H, c, mintc);
        qhlhopsize += H.size();
        //printf("%d %d", mtc, indH);
        for (int i = 0; i < H.size(); i++)
        {
            pathConcatenationC(L[s][H[i]].second, L[t][H[i]].second, C);
        }
    }
    //printf("%f\n", optw);
}

int main(int argc , char * argv[]){
    string s, sq;
    if (argc > 1)
        s = string(argv[1]);
    else
        s = string("NY");
    if (argc > 2)
        sq = string(argv[2]);
    else
        sq = string("q1");

    std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;

    t1=std::chrono::high_resolution_clock::now();
    load(s);
    t2=std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"Load Time "<<runT<<endl;


    string filename = string("../data/") + s + string("/") + sq;
    string answername = filename + string("QHLans");
    FILE *fp_query = fopen(filename.c_str(), "r");
    FILE *fp_ans = fopen(answername.c_str(), "w");
    t1=std::chrono::high_resolution_clock::now();
    int qs, qt;
    double qC;
    while (~fscanf(fp_query, "%d%d%lf", &qs, &qt, &qC)){
        QHL(qs, qt, qC);
        fprintf(fp_ans, "%f\n", optw);
    }
    fclose(fp_ans);
    t2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
    cout<<"QHL Query Time "<<runT<<endl;
    cout << "# of QHL Hoplinks " << qhlhopsize <<endl;
    cout << "# of QHL Path Concatenations " << qhltc <<endl;
    cout << "# hit Pruning Conditions " << hitpc << "in " << totvisitpcs <<endl;
}