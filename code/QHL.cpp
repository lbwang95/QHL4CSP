#include<bits/stdc++.h>
#include<unordered_map>
using namespace std;
typedef pair<double, double> DD;
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

void Cub_h(vector<DD> &P_, vector<DD> &P__, double &Ch){
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
unordered_map<int,vector<double>> Cub[MAX_V];//following ancarray positions -1
void pruningConditions(int cs, int v, bool build, vector<II> &H, double C){
    //printf("cs%d v%d\n", cs + 1, v + 1);
    if (build){
        if(Cub[cs].count(v) == 0){
            vector<II> _anc;
            curcs = cs;
            for (int i = 0; i < ancarray[cs].size() - 1; i++){//minus cs
                Cub[cs][v].push_back(0);
                _anc.push_back(II(ancarray[cs][i], i));
            }
            sort(_anc.begin(), _anc.end(), cmpp);
            int lena = _anc.size();
            for (int i = 0; i < lena; i++){
                for (int j = i + 1; j < lena; j++){
                    vector<DD> res;
                    int levu = _anc[i].first, levh = _anc[j].first;
                    int u = L[cs][levu].first, h = L[cs][levh].first;
                    concat(L[v][levu].second, (levu>levh)?L[u][levh].second:L[h][levu].second, res);
                    sort(res.begin(), res.end());
                    Cub_h(L[v][levh].second, res, Cub[cs][v][_anc[j].second]);
                    //printf("|%d %d %d %d %f|", levu, levh, u + 1, h + 1, Cub[cs][v][_anc[j].second]);
                }
            }
        }
    }
    if(Cub[cs].count(v) != 0){
        for (int i = 0; i < Cub[cs][v].size();i++){
            if(C>=Cub[cs][v][i])
                H.push_back(II(ancarray[cs][i], i));
        }
        //for (int i = 0; i < H.size();i++)
        //    printf("|%d %d|\n", H[i].first, H[i].second);
    }
}

void pathConcatenationC(vector<DD> &Psh,vector<DD> &Pht, double &C){
    int i = 0, j = Pht.size() - 1;
    while (i != Psh.size() && j != -1){
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

void QHL(int s, int t, double C){
    optw = DBL_MAX;
    if (s == t)
        return;
    s--;
    t--;
    vector<int> ancs,anct;
    int u1 = s, l = -1;
    while(u1!=rootpa){
        ancs.push_back(u1);
        u1 = T[u1].parent;
    }
    u1 = t;
    while(u1!=rootpa){
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
        int cs = ancs[i], ct = anct[j];
        //printf("cs%d ct%d\n", cs + 1, ct + 1);
        vector<II> H[4];
        pruningConditions(cs, s, 1, H[0], C);
        pruningConditions(cs, t, 1, H[1], C);
        pruningConditions(ct, s, 1, H[2], C);
        pruningConditions(ct, t, 1, H[3], C);
        //cost estimation
        int mtc = INT_MAX, indH;
        for (int i = 0; i < 4;i++){
            int timecost = 0;
            if(i<2){
                for (int j = 0; j < H[i].size(); j++){
                    int levh = ancarray[cs][H[i][j].second];
                    timecost += L[s][levh].second.size() + L[t][levh].second.size();
                }
            }
            else{
                for (int j = 0; j < H[i].size(); j++){
                    int levh = ancarray[ct][H[i][j].second];
                    timecost += L[s][levh].second.size() + L[t][levh].second.size();
                }
            }
            if(timecost<mtc){
                mtc = timecost;
                indH = i;
            }
        }
        npc += mtc;
        nhoplink += H[indH].size();
        //printf("%d %d", mtc, indH);
        if(indH<2){
            for (int i = 0; i < H[indH].size(); i++)
            {
                ind = ancarray[cs][H[indH][i].second];
                pathConcatenationC(L[s][ind].second, L[t][ind].second, C);
            }
        }
        else{
            for (int i = 0; i < H[indH].size(); i++)
            {
                ind = ancarray[ct][H[indH][i].second];
                pathConcatenationC(L[s][ind].second, L[t][ind].second, C);
            }
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
    cout << "Query Time " << runT << endl;
    cout << "# of path concatenations " << npc << endl;
    cout << "# of hoplinks " << nhoplink << endl;
}