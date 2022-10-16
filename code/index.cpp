#include<bits/stdc++.h>
#include<unordered_map>
using namespace std;
typedef pair<double, double> DD;
typedef pair<double, int> DI;
typedef pair<int, vector<DD>> IV;
typedef pair<int, vector<int>> IVI;
typedef pair<int, int> II;
typedef pair<int, double> ID;
const int MAX_V = 435676;
int N, M;
long long csp2hoptc, qhltc;
double optw = DBL_MAX;
int treeheight = 0, treewidth = 0, treeavgheight = 0;

typedef struct node{
    int level;
    int ran;
    int parent;
    vector<int> children;
    vector<IV> X;
};
node T[MAX_V];
vector<IV> L[MAX_V];
vector<int> anc[MAX_V];
int root = -1;
unordered_map<int,vector<DD> > adj[MAX_V];
unordered_map<int,int> adjo[MAX_V];
vector<int> order;
bool flag[MAX_V];
bool cmp(const IV &a, const IV &b){
    return T[a.first].ran > T[b.first].ran;
}

vector<int> ordergen;
int del[MAX_V];//deleted neighbors
double update(int v){
    return 1000 * adjo[v].size() + del[v];
}
typedef pair<II, int> III;
void genorder(string filename){
    priority_queue<II, vector<II>, greater<II> > degque;
    for (int i = 0; i < N; i++)
        degque.push(II(update(i), i));
    int iter = -1, totnewedge = 0;
    while(!degque.empty()){
        II ii = degque.top();
        degque.pop();
        int v = ii.second;
        if(flag[v])
            continue;
        double prio = update(v);
        if (prio > degque.top().first){
            degque.push(II(prio,v));
            continue;
        }
        iter += 1;
        //if(iter%100000==0){
        //    printf("----------%d %d--------\n", iter, totnewedge);
        //}
        flag[v] = 1;
        ordergen.push_back(v);
        T[v].ran = iter;
        unordered_map<int, int>::iterator it;
        vector<int> nei;
        for (it = adjo[v].begin(); it !=adjo[v].end(); it++)
            if(!flag[it->first])
                nei.push_back(it->first);
        int lenX = nei.size();
        for (int j = 0; j < lenX; j++){
            int u = nei[j];
            for (int k = j + 1; k < lenX; k++){
                int w = nei[k];
                if(adjo[u].count(w)==0){
                    adjo[u][w] = 1;
                    adjo[w][u] = 1;
                    totnewedge += 1;
                }
            }
            //adjo[u].erase(v);
            del[u]++;
        }
    }
    FILE *fp_order = fopen(filename.c_str(), "w");
    for (int i = 0; i < N;i++){
        fprintf(fp_order, "%d\n", T[i].ran);
    }
    fclose(fp_order);
}

void Sky(vector<DD> &P1,vector<DD> &P2, vector<DD> &res){
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
    //for (int i = 0; i < res.size();i++)
    //    printf("|%f %f|\n", res[i].first, res[i].second);
}

int descnt[MAX_V];
void treedec(){
    for (int i = 0; i < N; i++){
        int v = ordergen[i];
        //printf("|%d %d %d|\n", ii.first, adj[v].size(),v+1);
        if(i%100000==0)
            printf("%d\n", i);
        unordered_map<int, vector<DD>>::iterator it;  
        for (it = adj[v].begin(); it !=adj[v].end(); it++)
            T[v].X.push_back(IV(it->first, it->second));
        int lenX = T[v].X.size();
        //for (int j = 0; j < lenX;j++){
        //    printf("%d ", T[v].X[j].first+1);
        //}
        //cout << endl;
        for (int j = 0; j < lenX; j++){
            IV nu = T[v].X[j];
            int u = nu.first;
            for (int k = j + 1; k < lenX; k++){
                IV nw = T[v].X[k];
                int w = nw.first;
                vector<DD> emp;
                if(T[u].ran<T[w].ran){
                    if(adj[u].count(w)==0){
                        adj[u][w] = emp;
                        Sky(nu.second, nw.second, adj[u][w]);
                    }
                    else{
                        Sky(nu.second, nw.second, adj[u][w]);
                    }
                }
                else{
                    if(adj[w].count(u)==0){
                        adj[w][u] = emp;
                        Sky(nu.second, nw.second, adj[w][u]);
                    }
                    else{
                        Sky(nu.second, nw.second, adj[w][u]);
                    }
                }
            }
        }
    }
    for (int i = 0; i < ordergen.size();i++){
        int v = ordergen[i];
        sort(T[v].X.begin(), T[v].X.end(), cmp);
        int lenx = T[v].X.size();
        if(lenx!=0)
            T[v].parent = T[v].X[lenx - 1].first;
        else
            T[v].parent = MAX_V;
        //printf("%d %d\n", v + 1, T[v].parent + 1);
        vector<DD> tmpd;
        T[v].X.push_back(IV(v, tmpd));
        treewidth = max(treewidth, lenx + 1);
        if (T[v].parent == MAX_V){
            root = v;
            break;
        }
        T[T[v].parent].children.push_back(v);
        descnt[v]++;
        descnt[T[v].parent] += descnt[v];
    }
    cout << "Tree Width " << treewidth << endl;
}

queue<int> bfs;

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

vector<int> ancarray[MAX_V];
void generateLabel4v(int v){
    vector<int> anc;
    int u1 = v;
    //printf("v%d:", v+1);
    while(T[u1].parent!=MAX_V){
        anc.push_back(T[u1].parent);
        u1 = T[u1].parent;
    }
    int lenanc = anc.size();
    treeavgheight += lenanc;
    treeheight = max(treeheight, lenanc + 1);
    for (int i = 0; i < lenanc;i++){
        int u = anc[anc.size() - 1 - i];
        //printf("u%d ", u + 1);
        int lenx = T[v].X.size();
        vector<DD> res;
        for (int j = 0; j < lenx;j++){
            int w = T[v].X[j].first;
            if (w == v){
                continue;
            }
            //printf("w%d ", w + 1);
            if (T[w].ran <= T[u].ran){
                Sky(T[v].X[j].second, L[w][i].second, res);
                if(w==u)
                    ancarray[v].push_back(i);
            }
            else{//w>u, w has been in j-th ancarray //ancarray[v][j] not used because need sorting ancarray
                Sky(T[v].X[j].second, L[u][ancarray[v][j]].second, res);
            }
        }
        //for (int j = 0; j < res.size();j++)
        //    printf("res(%f,%f)", res[j].first, res[j].second);
        //vector<DD> tmp;
        //Sky(tmp, tmp, res, sres);
        //cout << endl;
        //for (int j = 0; j < sres.size();j++)
        //    printf("sres%d(%f,%f)", u+1, sres[j].first, sres[j].second);
        L[v].push_back(IV(u, res));
        //cout << endl;
    }
    vector<DD> tmpv;
    L[v].push_back(IV(v, tmpv));
    ancarray[v].push_back(anc.size());
    //for (int j = 0; j < ancarray[v].size();j++)
    //    printf("anc%d ", ancarray[v][j]);
    //cout << L[v].size() << endl;
}
void labeling(){
    bfs.push(root);
    int iter = 0;
    while(!bfs.empty()){
        int v= bfs.front();
        bfs.pop();
        //sort(T[v].X.begin(), T[v].X.end(), cmp);
        generateLabel4v(v);
        for (int i = 0; i < T[v].children.size();i++){
            bfs.push(T[v].children[i]);
        }
        if(iter%100000==0)
            printf("%d %d\n", iter, treeheight);
        iter += 1;
    }
}

queue<int> bfssave;
long long indexsize;
void save(string filename){
    filename += string("index");
    ofstream of;
    of.open(filename.c_str(), ios::binary);
    // FILE *fp_index=fopen("index.txt","w");
    // fprintf(fp_index, "%d ", N);
    of.write(reinterpret_cast<const char *>(&N), sizeof(int));
    bfssave.push(root);
    while(!bfssave.empty()){
        int v = bfssave.front();
        bfssave.pop();
        //printf("%d\n", v);
        int lenl = L[v].size(), nx = T[v].X.size();
        indexsize = indexsize + 4 + nx;
        //fprintf(fp_index, "%d %d %d %d%c", v, T[v].parent, nx, lenl,' ');
        of.write(reinterpret_cast<const char *>(&v), sizeof(int));
        of.write(reinterpret_cast<const char *>(&T[v].parent), sizeof(int));
        of.write(reinterpret_cast<const char *>(&nx), sizeof(int));
        of.write(reinterpret_cast<const char *>(&lenl), sizeof(int));
        for (int i = 0; i < nx; i++){
            //fprintf(fp_index, "%d%c", T[v].X[i].first, (i == nx - 1) ? ' ' : ' ');
            of.write(reinterpret_cast<const char *>(&T[v].X[i].first), sizeof(int));
        }
        for (int i = 0; i < lenl; i++){
            int lend = L[v][i].second.size();
            indexsize = indexsize + 2 + lend * 2;
            //fprintf(fp_index, "%d %d ", L[v][i].first, lend);
            of.write(reinterpret_cast<const char *>(&L[v][i].first), sizeof(int));
            of.write(reinterpret_cast<const char *>(&lend), sizeof(int));
            for (int j = 0; j < lend;j++){
                //fprintf(fp_index, "%d %d ", L[v][i].second[j].first*10, L[v][i].second[j].second*10);
                int p = L[v][i].second[j].first, q = L[v][i].second[j].second;
                of.write(reinterpret_cast<const char *>(&p), sizeof(int));
                of.write(reinterpret_cast<const char *>(&q), sizeof(int));
            }
        }
        //fprintf(fp_index, "\n");
        //sort(T[v].X.begin(), T[v].X.end(), cmp);
        for (int i = 0; i < T[v].children.size();i++){
            bfssave.push(T[v].children[i]);
        }
    }
    //fclose(fp_index);
    of.close();
}

int curcs;
bool cmpp(const II &a,const II &b){
    return L[curcs][a.first].second[0].first < L[curcs][b.first].second[0].first;
}

void Cub_h(int v,int u,int h,vector<DD> &P_, vector<DD> &P__, double &Ch){
    /*
    for (int k = 0; k < P_.size();k++)
        printf("*%f %f*", P_[k].first, P_[k].second);
    for (int k = 0; k < P__.size();k++)
        printf("#%f %f#", P__[k].first, P__[k].second);*/
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
unordered_map<int, double> relationC[MAX_V];
unordered_map<int,vector<DI>> Cub[MAX_V];//following ancarray positions -1
int totpcs = 0, totvisitpcs = 0;
long long qhlindexsize;
void pruningConditions(int cs, int v){
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
    /*
    for (int i = 0; i < lena; i++){
        for (int j = i + 1; j < lena; j++){
            vector<DD> res;
            int levu = _anc[i].first, levh = _anc[j].first;
            int u = L[cs][levu].first, h = L[cs][levh].first;
            if(relationC[v].count(levu*treeheight+levh)!=0){
                Cub[v][cs][_anc[j].second].first=max(Cub[v][cs][_anc[j].second].first,relationC[v][levu*treeheight+levh]);
            }
            else{
                Sky(L[v][levu].second, (levu>levh)?L[u][levh].second:L[h][levu].second, res);
                Cub_h(v,u,h,L[v][levh].second, res, Cub[v][cs][_anc[j].second].first);
                relationC[v][levu * treeheight + levh] = Cub[v][cs][_anc[j].second].first;
                //printf("|%d %d %d %d %f|", levu, levh, u + 1, h + 1, Cub[v][cs][_anc[j].second]);
            }
        }
    }*/
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

void QHLindex(string prefix){
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

void pathConcatenationC(vector<DD> &Psh,vector<DD> &Pht, double &C){
    int i = 0, j = Pht.size() - 1;
    while (i != Psh.size() && j != -1){
        qhltc++;
        if (Psh[i].first + Pht[j].first <= C){
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
int qhlhopsize, csp2hopsize;
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

void concatQuery(vector<DD> &P1,vector<DD> &P2, double C){
    for (int i = 0; i < P1.size();i++)
        for (int j = 0; j < P2.size();j++){
            csp2hoptc++;
            //printf("(%f,%f-%f,%f)", P1[i].first, P1[i].second, P2[j].first, P2[j].second);
            if (P1[i].first + P2[j].first <= C){
                optw=min(optw, P1[i].second + P2[j].second);
            }
        }
}

void CSP2Hop(int s, int t, double C){
    optw = DBL_MAX;
    if (s == t)
        return;
    s--;
    t--;
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
            i--;
            j--;
            k++;
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
        //printf("lca%d\n", l + 1);
        csp2hopsize += ancarray[l].size();
        for (int i = 0; i < ancarray[l].size();i++){
            ind = ancarray[l][i];
            //printf("|%d %d %d|\n", ind, L[s][ind].second.size(), L[t][ind].second.size());
            concatQuery(L[s][ind].second, L[t][ind].second, C);
        }
    }
    //printf("%f\n", optw);
}
struct edge{
    int from, to;
    double cost, weight;
    edge(int a,int b,double c,double d){
        from = a, to = b, cost = c, weight = d;
    }
};
vector<edge> alledges;

int main(int argc , char * argv[]){
    string s, sq;
    FILE *fp_query, *fp_networkw, *fp_networkc;
    if (argc > 1)
        s = string(argv[1]);
    else
        s = string("NY");
    if (argc > 2)
        sq = string(argv[2]);
    else
        sq = string("pool1");
    string prefix = string("../data/") + s + string("/");
    string s1 = prefix + string("USA-road-t.") + s + (".gr");
    string s2 = prefix + string("USA-road-d.") + s + (".gr");
    //test or not
    if(0){
        sq = string("test");
        s1 = string("../data/") + s + string("/t.") + s + (".gr");
        s2 = string("../data/") + s + string("/d.") + s + (".gr");
    }
    fp_networkw = fopen(s1.c_str(), "r");
    fp_networkc = fopen(s2.c_str(), "r");
    char ch, buffer[100];
    int u, v;
    double w, c;
    //travel time
    for (int i = 0; i < 4; i++)
        fgets(buffer, 90, fp_networkw);
    for (int i = 0; i < 4; i++)
        fgetc(fp_networkw);
    fscanf(fp_networkw, "%d%d", &N, &M);
    for (int i = 0; i < 3; i++)
        fgets(buffer, 90, fp_networkw);
    //distance
    for (int i = 0; i < 4; i++)
        fgets(buffer, 90, fp_networkc);
    for (int i = 0; i < 4; i++)
        fgetc(fp_networkc);
    fscanf(fp_networkc, "%d%d", &N, &M);
    for (int i = 0; i < 3; i++)
        fgets(buffer, 90, fp_networkc);
    for (int i = 0; i < M; i++) {
        fscanf(fp_networkw, "%c%d%d%lf", &ch, &u, &v, &w);
        fgets(buffer, 90, fp_networkw);
        fscanf(fp_networkc, "%c%d%d%lf", &ch, &u, &v, &c);
        fgets(buffer, 90, fp_networkc);
        u--;
        v--;        
        //printf("%d %d\n", u, v);
        if(i%2==0){
            adjo[u][v]=1;
            adjo[v][u]=1;
            alledges.push_back(edge(u, v, c, w));
        }
    }

    std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;
    t1=std::chrono::high_resolution_clock::now();
    string ordername = string("../data/") + s + string("/") + string("order.txt");
    //genorder(ordername);
    ordergen.assign(N, -1);
    FILE *fpord = fopen(ordername.c_str(), "r");
    for (int i = 0; i < N; i++){
        fscanf(fpord, "%d", &T[i].ran);
        ordergen[T[i].ran] = i;
    }
    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"Order Generation Time "<<runT<<endl;
    
    //distribute edges
    for (int i = 0; i < alledges.size();i++){
        int f = alledges[i].from, t = alledges[i].to;
        if(T[f].ran<T[t].ran)
            adj[f][t].push_back(DD(alledges[i].cost, alledges[i].weight));
        else
            adj[t][f].push_back(DD(alledges[i].cost, alledges[i].weight));
    }

    t1=std::chrono::high_resolution_clock::now();
    treedec();
    t2=std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"Tree Decomposition Time "<<runT<<endl;
    //return 0;

    t1=std::chrono::high_resolution_clock::now();
    labeling();
    t2=std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"Labeling Time "<<runT<<endl;
    cout<<"Tree Height "<<treeheight<<endl;
    cout<<"Tree Avg. Height "<<(double)treeavgheight/N<<endl;

    t1=std::chrono::high_resolution_clock::now();
    QHLindex(prefix);
    t2=std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"QHL Index Time "<<runT<<endl;
    cout << "QHL Index Size " << (double)qhlindexsize * 8 / 1000000 << "MB" << endl;

    
    t1=std::chrono::high_resolution_clock::now();
    //FILE *fp = fopen("index", "w");
    //save(string("../data/")+s+string("/"));
    t2=std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	//cout<<"Saving Time "<<runT<<endl;
    cout << "Index Size " << (double)indexsize * 4 / 1000000 << "MB" << endl;
    freopen((prefix + string("Results")).c_str(), "w", stdout);
    for (int i = 0; i < 5;i++){
        totvisitpcs = hitpc = qhlhopsize = csp2hopsize = qhltc = csp2hoptc = 0;
        string s3 = string("../data/") + s + string("/") + string("q") + to_string(i+1);
        fp_query = fopen(s3.c_str(), "r");
        vector<pair<II, double>> queryset;
        int qs, qt;
        double qC;
        while (~fscanf(fp_query, "%d%d%lf", &qs, &qt, &qC)){
            queryset.push_back(make_pair(II(qs, qt), qC));
        }
        FILE *fp_ans1 = fopen((s3+string("CSP2Hopans")).c_str(), "w");
        t1=std::chrono::high_resolution_clock::now();
        for(int i=0;i<queryset.size();i++){
            CSP2Hop(queryset[i].first.first, queryset[i].first.second, queryset[i].second);
            fprintf(fp_ans1, "%f\n", optw);
        }
        fclose(fp_ans1);
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT= time_span.count();
        cout<<"CSP-2Hop Query Time "<<runT<<endl;
        cout << "# of CSP-2Hop Hoplinks " << csp2hopsize <<endl;
        cout << "# of CSP-2Hop Path Concatenations " << csp2hoptc <<endl;

        FILE *fp_ans2 = fopen((s3+string("QHLans")).c_str(), "w");
        t1=std::chrono::high_resolution_clock::now();
        for(int i=0;i<queryset.size();i++){
            QHL(queryset[i].first.first, queryset[i].first.second, queryset[i].second);
            fprintf(fp_ans2, "%f\n", optw);
        }
        fclose(fp_ans2);
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT= time_span.count();
        cout<<"QHL Query Time "<<runT<<endl;
        cout << "# of QHL Hoplinks " << qhlhopsize <<endl;
        cout << "# of QHL Path Concatenations " << qhltc <<endl;
        cout << "# hit Pruning Conditions " << hitpc << "in " << totvisitpcs <<endl;
    }
    for (int i = 0; i < 5;i++){
        totvisitpcs = hitpc = qhlhopsize = csp2hopsize = qhltc = csp2hoptc = 0;
        string s3 = string("../data/") + s + string("/") + string("r") + to_string(i+1);
        fp_query = fopen(s3.c_str(), "r");
        vector<pair<II, double>> queryset;
        int qs, qt;
        double qC;
        while (~fscanf(fp_query, "%d%d%lf", &qs, &qt, &qC)){
            queryset.push_back(make_pair(II(qs, qt), qC));
        }
        FILE *fp_ans1 = fopen((s3+string("CSP2Hopans")).c_str(), "w");
        t1=std::chrono::high_resolution_clock::now();
        for(int i=0;i<queryset.size();i++){
            CSP2Hop(queryset[i].first.first, queryset[i].first.second, queryset[i].second);
            fprintf(fp_ans1, "%f\n", optw);
        }
        fclose(fp_ans1);
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT= time_span.count();
        cout<<"CSP-2Hop Query Time "<<runT<<endl;
        cout << "# of CSP-2Hop Hoplinks " << csp2hopsize <<endl;
        cout << "# of CSP-2Hop Path Concatenations " << csp2hoptc <<endl;

        FILE *fp_ans2 = fopen((s3+string("QHLans")).c_str(), "w");
        t1=std::chrono::high_resolution_clock::now();
        for(int i=0;i<queryset.size();i++){
            QHL(queryset[i].first.first, queryset[i].first.second, queryset[i].second);
            fprintf(fp_ans2, "%f\n", optw);
        }
        fclose(fp_ans2);
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT= time_span.count();
        cout<<"QHL Query Time "<<runT<<endl;
        cout << "# of QHL Hoplinks " << qhlhopsize <<endl;
        cout << "# of QHL Path Concatenations " << qhltc <<endl;
        cout << "# hit Pruning Conditions " << hitpc << "in " << totvisitpcs <<endl;
    }
}