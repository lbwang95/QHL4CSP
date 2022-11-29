#include <cstdio>
#include "utility1.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <queue>
#include <cassert>
#include <algorithm>
#include <map>
#include <time.h>
#include <cmath>
#include <set>
#include <queue>
#include <climits>
#include <cfloat>
#include <stack>
#include <queue>
#include <unordered_map>
#include <omp.h>
#include <iomanip>


bool pair_comp_func_long_bigfirst( pair<int, long> p1, pair<int, long> p2){
    return p1.second > p2.second;
}

bool lambda_tree_sorter( pair<double, pair<int, int> > p1, pair<double, pair<int, int> > p2){
    return p1.first > p2.first;
}

bool single_dimension_sorter( pair<int, pair<int, int> > p1, pair<int, pair<int, int> > p2){
    if(p1.second.first > p2.second.first) return true;
    else if (p1.second.first < p2.second.first) return false;
    else if( p1.second.second > p2.second.second) return true;
    else return false;
}



using namespace std;
const int VectorDefaultSize=1;
const int TOPNUM = 1;
const int BOUNDARY_QUERY_NUM=1000000;
const int QUERY_NUM=1000000;

const int PARTITION_SIZE = 2048;

//#define SQRT_APPROX 1.0488
#define BOUNDARY_GRAPH_APPROX 1.0
#define SQRT_APPROX 1.0
#define APPROX 1.0

struct query{
    int u;
    int v;
    int cost_limit;
};

struct edge_list{
    int u;
    int v;
    int cost;
    int weight;
    edge_list(int uu, int vv, int ccost, int wweight){
        u = uu;
        v = vv;
        cost = ccost;
        weight = wweight;
    }

    bool operator ==( const edge_list& e) const{
        return u== e.u && v == e.v && cost == e.cost && weight == e.weight;
    }

    edge_list(){}
};

struct link_list_edge{
    int v;
    int weight;
    int cost;
    //int dominated_cost;
    int dominate_paths_minimum_weight;
    bool operator<(const link_list_edge &tri) const{
        if(cost > tri.cost) return true;
        else if(cost < tri.cost ) return false;
        else if(weight > tri.weight ) return true;
        else if(weight < tri.weight) return false;
        else if ( v > tri.v ) return true;
        else return false;
    }
    
    link_list_edge(int vv, int wweight, int ccost){
        v = vv;
        weight = wweight;
        cost = ccost;
        dominate_paths_minimum_weight = wweight;
        //dominated_cost = ccost;
    }
    link_list_edge(){}
};

bool paretor_path_compartor( pair<pareto_pair<int, int>, pair<int,int> > p1, pair<pareto_pair<int, int>, pair<int, int> > p2){
    if(p1.first.first > p2.first.first) return true;
    else if( p1.first.first < p2.first.first ) return false;
    else if(p1.first.second > p2.first.second ) return true;
    else if(p1.first.second < p2.first.second) return false;
    else if( p1.second > p2.second ) return true;
    else return false;
}


template <typename __T1>
struct triple{
    int v;
    int weight;
    __T1 cost;
    int pivot_level;
    int left_child_first_pointer;
    int left_child_second_pointer;
    int right_child_first_pointer;
    int right_child_second_pointer;
    int dominated_weight;
    int dominated_cost;
    bool operator<(const triple &tri) const{
        if(cost > tri.cost) return true;
        else if(cost < tri.cost ) return false;
        else if(weight > tri.weight ) return true;
        else if(weight < tri.weight) return false;
        else if ( v > tri.v ) return true;
        else return false;
    }

    bool operator==(const triple &tri)const{
        return v == tri.v && weight == tri.weight && 
            cost == tri.cost;
    }

    triple(int vv, int wweight, __T1 ccost){
        v = vv;
        weight = wweight;
        cost = ccost;
        pivot_level = -2;
        left_child_first_pointer = -5;
        left_child_second_pointer = -5;
        right_child_first_pointer = -5;
        right_child_second_pointer = -5;
        dominated_weight = wweight;
        dominated_cost = ccost;
    }
    triple(){}
};

bool triple_small_first(triple<int> t1, triple<int> t2){
    if(t1.cost < t2.cost ) return true;
    else if(t1.cost > t2.cost) return false;
    else if(t1.weight < t2.weight) return true;
    else  return false;
}


template <typename __T1>
struct label_cache{
    int v;
    int min_cost;
    int start_pos;
    int end_pos;
};

template <typename __T1>
struct label_structure{
    iVector<label_cache<__T1> > label_nodes;
    iVector<triple<__T1> > label_content;
};


struct path_hash_keys{
    int u;
    int v;
    int weight;
    int cost;
    int direction;
    path_hash_keys(int uu, int vv, int costt, int weightt, int directionn){
         u = uu;
         v = vv;
         cost = costt;
         weight = weightt;
         direction = directionn;
    }
    path_hash_keys(){}

    bool operator<(const path_hash_keys &hk) const{
        if(u < hk.u) return true;
        else if( u > hk.u) return false;
        else if( v < hk.v ) return true;
        else if( v > hk.v) return false;
        else if( cost < hk.cost) return true;
        else if( cost > hk.cost) return false;
        else if(weight < hk.weight) return true;
        else if(weight > hk.weight) return false;
        else if(direction < hk.direction) return true;
        else return false;
    }
    bool operator==(const path_hash_keys &hk)const{
        return u==hk.u && v == hk.v && cost == hk.cost && weight == hk.weight && direction == hk.direction;
    }

    bool operator=(const path_hash_keys & hk){
        u = hk.u;
        v = hk.v;
        cost = hk.cost;
        weight = hk.weight;
        direction = hk.direction;
    }
};


struct CSP{
    int N;
    int M;
    int BG_N;
    int BG_M;
    int* level;
    int* v2level;
    bool* pruned;
    int* partition;
    int approximate_reduced_labels;
    int total_labels;
    int max_cost;
    int avg_cost;
    int* cost_min_query;
    int* weight_max_query;
    int* weight_min_query;
    int* cost_min_query_2;
    int* weight_max_query_2;
    int* weight_min_query_2;
    int num_of_partitions;
    double generate_boundary_time;
    double ordering_time;
    double index_construction_time;
    iMap<int> vid2bid;
    iMap<int> bid2vid;
    //query* boundary_query_set;
    query* query_set;
    pair<int, int>* last_pruned;
    iVector<vector<link_list_edge > > original_links[2];
    iVector<vector<link_list_edge > > links[2];
   
    iVector<label_structure<int> > opt_labels[2];

    iVector<vector<vector<triple<int> > > > labels[2];
    iVector<vector<int> > label_nodes[2]; 
    iVector<iVector<int> > partition_boundary;
    std::priority_queue<traversal_tuple<int> > pq;
    iHeap<int> bidij_1d[2];
    std::priority_queue<traversal_tuple<int> >* individual_pq;
    MMap<vector<pareto_pair<int, int> > > mark[2];
    MMap<map<struct path_hash_keys, int> > pre_paths[2]; 
   
    iVector<map<int, int> > tree_father;
    iVector<map<int, vector<int> > > tree_children;
    iVector<map<int, int> > tree_coverage;

    MMap<vector<triple<int> > > label_mark;
    MMap<int> start_prunning;
    iMap<int> pq_mark;
    iMap<int> pq_mark_back;
    MMap<int> pq_mark_max_weight;
    MMap<int> pq_mark_min_weight;
    Traversal_Heap tra_pq;
    MMap<traversal_tuple<int> > tra_estimate;

    void init_for_G(){

        original_links[0].re_allocate(N);
        original_links[0].m_num = N;
        original_links[1].re_allocate(N);
        original_links[1].m_num = N;

        vector<triple<int> > nil_lbl;
        label_mark.set_nil(nil_lbl);
        label_mark.initialize(N);
        label_mark.clean();

        generate_boundary_time = 0.0;
        ordering_time = 0.0;
        index_construction_time = 0.0;

        //boundary_query_set = new query[BOUNDARY_QUERY_NUM];
        bidij_1d[0].initialize(N);
        bidij_1d[1].initialize(N);


        tra_pq.initialize(N);
        traversal_tuple<int> nil_tp(-1,-1,-1);
        tra_estimate.set_nil(nil_tp);
        tra_estimate.initialize(N);
        tra_estimate.clean();
        vid2bid.initialize(N); 
        bid2vid.initialize(N); 
        
        query_set = new query[QUERY_NUM];
        level = new int[N];
        v2level = new int[N];
        last_pruned = new pair<int, int>[N];
        for( int x=0; x<N; x++){
            last_pruned[x].first = -1;
            last_pruned[x].second = -1;
        }
        pruned = new bool[N];
        partition = new int[N];

        max_cost = -1;
        avg_cost = 0;
        mark[0].initialize(N);
        mark[1].initialize(N);
        mark[0].clean();
        mark[1].clean();
        pq_mark.initialize(N);
        pq_mark.clean();
        pq_mark_back.initialize(N);
        pq_mark_back.clean();
        int nil_int = -10;
        pq_mark_max_weight.set_nil(nil_int);
        pq_mark_max_weight.initialize(N);
        pq_mark_min_weight.set_nil(nil_int);
        pq_mark_min_weight.initialize(N);
        start_prunning.set_nil(nil_int);
        start_prunning.initialize(N);

        approximate_reduced_labels = 0;
        total_labels=0;
    }

    void init_for_BG(){
        labels[0].re_allocate(BG_N);
        labels[0].m_num = BG_N;
        labels[1].re_allocate(BG_N);
        labels[1].m_num = BG_N;
        label_nodes[0].re_allocate(BG_N);
        label_nodes[1].re_allocate(BG_N);
        label_nodes[0].m_num = BG_N;
        label_nodes[1].m_num = BG_N;
        opt_labels[0].re_allocate(BG_N);
        opt_labels[1].re_allocate(BG_N);
        opt_labels[0].m_num = BG_N;
        opt_labels[1].m_num = BG_N;
        for(int i=0; i<BG_N; i++){
            opt_labels[0][i].label_content.clean();
            opt_labels[0][i].label_nodes.clean();
            opt_labels[1][i].label_content.clean();
            opt_labels[1][i].label_nodes.clean();
        }
        cost_min_query = new int[BG_N];
        weight_max_query = new int[BG_N];
        weight_min_query = new int[BG_N];
        cost_min_query_2 = new int[BG_N];
        weight_max_query_2 = new int[BG_N];
        weight_min_query_2 = new int[BG_N];
        links[0].re_allocate(BG_N);
        links[0].m_num = BG_N;
        links[1].re_allocate(BG_N);
        links[1].m_num = BG_N;
        individual_pq = new priority_queue<traversal_tuple<int> >[BG_N];
        for(int x=0; x<BG_N; x++){
            assert(individual_pq[x].empty());
        }

        map<struct path_hash_keys, int> nil_path;
        pre_paths[0].set_nil(nil_path);
        pre_paths[1].set_nil(nil_path);
        pre_paths[0].initialize(BG_N);
        pre_paths[1].initialize(BG_N);
        pre_paths[0].clean();
        pre_paths[1].clean();
    }

    void load_directed_graph(string graph_name){
        clock_t load_graph_start = clock();
        FILE* file = fopen(graph_name.c_str(), "r");
        int ret_val = fscanf(file, "%d", &N);
        if(ret_val <0){ cout<<"read input graph file error"<<endl; exit(0);}
        cout<<N<<" nodes"<<endl;
        BG_N = N;
        links[0].re_allocate(BG_N);
        links[0].m_num = BG_N;
        links[1].re_allocate(BG_N);
        links[1].m_num = BG_N;
        init_for_G();
        init_for_BG();
        int u, v, weight, cost;
        int lines =0;
        double tmp_cost = 0;
        M=0;
        for(; fscanf(file, "%d:", &u)==1; lines++){
            for(; fscanf(file, "%d", &v)==1;){
                if(v==-1) break;
                int ret_val = fscanf(file, "%d %d", &cost, &weight);
                if(ret_val <0){cout<<"input file error"<<endl; exit(0);}
                link_list_edge forward_tri(v, weight, cost);  
                link_list_edge backward_tri(u, weight, cost);
                if(cost > max_cost) max_cost = cost;
                links[0][u].push_back(forward_tri);
                links[1][v].push_back(backward_tri);
                tmp_cost=  M*1.0/(M+1) * tmp_cost +cost*1.0/(M+1);
                M++; 
            }
        }
        BG_M = M;
        avg_cost = tmp_cost;
        cout<<"Avg cost:"<<avg_cost<<endl;
        cout<<"number of edges: "<<M<<endl;
        mark[0].clean();
        fclose(file);
        clock_t load_graph_stop = clock();
        printf( "Load Graph Time: %0.2lf\n",(double)(load_graph_stop-load_graph_start)/CLOCKS_PER_SEC );
    }

    void load_directed_graph(string graph_name, string partition_name, string cut_name){
        clock_t load_graph_start = clock();
        FILE* file = fopen(graph_name.c_str(), "r");
        FILE* cut_file = fopen(cut_name.c_str(), "r");
        ifstream partition_file(partition_name);
        int ret_val = fscanf(file, "%d", &N);
        if(ret_val <0){ cout<<"read input graph file error"<<endl; exit(0);}
        cout<<N<<" nodes"<<endl;
        init_for_G();
        int u, v, weight, cost;
        int lines =0;
        double tmp_cost=0;
        M=0;
        for(; fscanf(file, "%d:", &u)==1; lines++){
            for(; fscanf(file, "%d", &v)==1;){
                if(v==-1) break;
                int ret_val = fscanf(file, "%d %d", &cost, &weight);
                if(ret_val <0){cout<<"input file error"<<endl; exit(0);}
                //cout<<u<<"->"<<v<<" "<<weight<<" "<<cost<<endl;
                link_list_edge forward_tri(v, weight, cost);  
                link_list_edge backward_tri(u, weight, cost);
                if(cost > max_cost) max_cost = cost;
                original_links[0][u].push_back(forward_tri);
                original_links[1][v].push_back(backward_tri);
                tmp_cost=  M*1.0/(M+1) * tmp_cost +cost*1.0/(M+1);
                M++;
            }
        }
        avg_cost= (int)tmp_cost;
        cout<<"Avg cost:"<<avg_cost<<endl;
        cout<<"number of edges: "<<M<<endl;
      
        clock_t load_graph_stop = clock();
        printf( "Load Graph Time: %0.2lf\n",(double)(load_graph_stop-load_graph_start)/CLOCKS_PER_SEC );
        

        for(int xx=0; xx<N; xx++){
            pruned[xx]=1;
            partition[xx] =-1;
        }


        cout<<"reading partition file"<<endl;
        string buf;
        string delimt("============================");
        int current_part=-1;
        while(getline(partition_file, buf)){
            buf.erase(buf.find_last_not_of("\n\r\t")+1); 
            if(buf.compare(delimt) ==0 ){
                getline(partition_file, buf);
                buf.erase(buf.find_last_not_of(": \n\r\t")+1); 
                buf = buf.substr(10);
                int partition_num = atoi(buf.c_str());
                //cout<<partition_num<<endl;
                current_part = partition_num;
            }else{
                stringstream ss(buf);
                int node_id;
                double x_axis, y_axis;
                ss>>node_id>>x_axis>>y_axis;
                node_id--;
                partition[node_id] = current_part;
            }
        }

        for(int x=0; x<N; x++){
            if(partition[x]<0) assert(0);
        }

        cout<<"finished reading partition file"<<endl;
        cout<<"number of partitions:"<<(current_part+1)<<endl;

        partition_boundary.re_allocate(current_part+1);
        partition_boundary.m_num = current_part+1;

        iVector<pair<int,int> > cut_vector;
        int cut_start_v, cut_end_v;
        for(lines =0; fscanf(cut_file, "%d %d", &cut_start_v, &cut_end_v)==2; lines++){
            
            pair<int, int> pr;
            if(cut_start_v < cut_end_v){
                pr.first = cut_start_v;
                pr.second = cut_end_v;
            }else{
                pr.first = cut_end_v;
                pr.second = cut_start_v;
            }
            cut_vector.push_back(pr);
        }

        cut_vector.unique();

        for(int x=0; x<cut_vector.size(); x++){
            cut_start_v = cut_vector[x].first;
            cut_end_v = cut_vector[x].second;
            pruned[cut_start_v]=0;
            //pruned[cut_end_v]=0;
            partition_boundary[partition[cut_start_v]].push_back(cut_start_v);
            partition_boundary[partition[cut_end_v]].push_back(cut_start_v);
        }

        for(int x=0; x<partition_boundary.m_num; x++){
            //cout<<"The "<<x<<"'th partition: "<<partition_boundary[x].size()<<endl;
            partition_boundary[x].unique();
        }



        vid2bid.initialize(N);
        int id=0;
        for(int i=0; i<N; i++){
            if(!pruned[i]){
                vid2bid.insert(i, id);
                id++;
            }
        }
        BG_N = id;

        bid2vid.initialize(BG_N);

        for(int i=0; i<vid2bid.occur.size(); i++){
            int x = vid2bid.occur[i];
            if(vid2bid.exist(x)){
                int bid = vid2bid[x];
                bid2vid.insert(bid, x);
            }
        }


        init_for_BG();


        string contrat_file_name = graph_name;
        contrat_file_name = "bg."+contrat_file_name +"a"+to_string(APPROX) ;

        FILE* contrat_file = fopen( contrat_file_name.c_str(), "w");
        
        fprintf(contrat_file, "%d\n", BG_N);


        int pruned_count=0;

        clock_t generate_boundary_start = clock();
        //int pruned_shortcut_count=0;
        for(int i=0; i<N; i++){
            if((i+1)%500000 == 0){
                mark[0].initialize(N);
                //bidij_1d[0].DeepClean();
            }
            //if(i%10000==0) cout<<i<<endl;
            if(pruned[i]){pruned_count++; continue;}
            mark[0].clean();
            //cout<<i<<endl;
            create_boundary_edges(i, contrat_file);
        }
        clock_t generate_boundary_end = clock();

         generate_boundary_time = (double)(generate_boundary_end-generate_boundary_start)/CLOCKS_PER_SEC ;
        cout<<"overlay graph: N="<<BG_N<<",  M="<<BG_M<<endl;
        cout<<"generate overlay time:"<<generate_boundary_time<<endl;
        mark[0].initialize(N);
        fclose(contrat_file);

        fclose(file);
        fclose(cut_file);
    }


    void in_boundary_single_dimension_by_cost(int u, int dir, int cost_limit){
       pq_mark.clean();
       pq_mark_max_weight.clean();
       pq_mark.insert(u, 0);
       pq_mark_max_weight.insert(u, 0);
       bidij_1d[dir].clean();
       bidij_1d[dir].insert(u, 0);
       while(!bidij_1d[dir].empty()){
            int forward_v  = bidij_1d[dir].head();
            int cost = pq_mark[forward_v];
            int weight = pq_mark_max_weight[forward_v];
            bidij_1d[dir].pop();
            if(!pruned[forward_v] && forward_v != u){
                continue;
            }
            for(int i=0; i<original_links[dir][forward_v].size(); i++){
                int expand_v = original_links[dir][forward_v][i].v;
                int expand_cost = original_links[dir][forward_v][i].cost;
                int expand_weight = original_links[dir][forward_v][i].weight;
                int total_cost = cost + expand_cost;
                int total_weight = weight + expand_weight;
                if(total_cost > cost_limit) continue;
                
                if(pq_mark.exist(expand_v)){
                    if(total_cost < pq_mark[expand_v]){
                        pq_mark.insert( expand_v, total_cost);
                        pq_mark_max_weight.insert( expand_v, total_weight);
                        bidij_1d[dir].insert(expand_v, total_cost);
                    }else if (total_cost == pq_mark[expand_v]){
                        if(pq_mark_max_weight[expand_v] > total_weight){
                            pq_mark_max_weight[expand_v] = total_weight;
                        }
                    }
                }else{
                    pq_mark.insert( expand_v, total_cost );
                    pq_mark_max_weight.insert( expand_v, total_weight);
                    bidij_1d[dir].insert(expand_v, total_cost);
                }
            }
       }
    }

    void in_boundary_single_dimension_by_weight(int u, int dir, int cost_limit){
        pq_mark.clean();
        pq_mark.insert(u, 0);
        pq_mark_min_weight.clean();
        pq_mark_min_weight.insert(u, 0);
        bidij_1d[dir].clean();
        bidij_1d[dir].insert(u, 0);
        while(!bidij_1d[dir].empty()){
            int forward_v  = bidij_1d[dir].head();
            int weight = pq_mark_min_weight[forward_v];
            int cost = pq_mark[forward_v];
            bidij_1d[dir].pop();
            if(!pruned[forward_v] && forward_v != u){
                continue;
            }
            for(int i=0; i<original_links[dir][forward_v].size(); i++){
                int expand_v = original_links[dir][forward_v][i].v;
                int expand_weight = original_links[dir][forward_v][i].weight;
                int expand_cost = original_links[dir][forward_v][i].cost;
                int total_weight = weight + expand_weight;
                int total_cost = cost+ expand_cost;
                if(total_cost > cost_limit) continue;
                if(pq_mark_min_weight.exist(expand_v)){
                    if(total_weight < pq_mark_min_weight[expand_v]){
                        pq_mark_min_weight.insert( expand_v, total_weight);
                        bidij_1d[dir].insert(expand_v, total_weight);
                        pq_mark.insert(expand_v, total_cost);
                    }
                }else{
                    pq_mark_min_weight.insert( expand_v, total_weight );
                    bidij_1d[dir].insert(expand_v, total_weight);
                    pq_mark.insert(expand_v, total_cost);
                }
            }
        }
    }

    void create_boundary_edges(int i, FILE* contrat_file){
        start_prunning.clean();

        in_boundary_single_dimension_by_weight(i, 0, INT_MAX);
        in_boundary_single_dimension_by_cost(i, 0, INT_MAX);

        for(int k=0; k<pq_mark_min_weight.occur.size(); k++){
            int reachable_v = pq_mark_min_weight.occur[k];
            int min_w = pq_mark_min_weight[reachable_v];
            int max_w = pq_mark_max_weight[reachable_v];
            if(min_w>0 && max_w >0){
                double val = max_w*1.0/min_w;
                int number_of_stored_paths = 1+(int)(log(val)/log(APPROX));
                start_prunning.insert(reachable_v, number_of_stored_paths);
            }
        }

        vector<traversal_tuple<int> > t_list;
        traversal_tuple<int> tt(i, 0, 0);
        tt.dominate_paths_minimum_weight = 0;
        pq.push(tt);
        while(!pq.empty()){
            traversal_tuple<int> top_t = pq.top();
            pq.pop();
            int to_v = top_t.v;
            int to_c = top_t.cost;
            int to_w = top_t.weight;
            int dominate_paths_minimum_weight = top_t.dominate_paths_minimum_weight;
            pareto_pair<int, int> pp(to_c, to_w);
            if(mark[0].exist(to_v)){
                int l_sz = mark[0][to_v].size();

                if(mark[0][to_v][l_sz-1].first <= to_c 
                   && mark[0][to_v][l_sz-1].dominate_paths_minimum_weight 
                   <=dominate_paths_minimum_weight){
                    continue;
                }else{
                    if(l_sz+1 > start_prunning[to_v]){
                        pp.dominate_paths_minimum_weight = max( pq_mark_min_weight[to_v],  int (ceil(to_w/BOUNDARY_GRAPH_APPROX) ) );
                    }else{
                        pp.dominate_paths_minimum_weight = dominate_paths_minimum_weight;
                    }
                    mark[0][to_v].push_back(pp);
                }
                if(!pruned[to_v] && to_v != i) continue;
            }else{
                if(start_prunning[to_v] > 1){
                    pp.dominate_paths_minimum_weight = max( pq_mark_min_weight[to_v],  int (ceil(to_w/BOUNDARY_GRAPH_APPROX) ) );
                }else{
                    pp.dominate_paths_minimum_weight = dominate_paths_minimum_weight;
                }
                vector<pareto_pair<int, int> > p_list;
                p_list.push_back(pp);
                mark[0].insert(to_v, p_list);
                if(!pruned[to_v] && to_v != i) continue;
            }
            dominate_paths_minimum_weight = pp.dominate_paths_minimum_weight;

            for(int x=0; x<original_links[0][to_v].size(); x++){
                int e_v = original_links[0][to_v][x].v; 
                int e_c = original_links[0][to_v][x].cost; 
                int e_w = original_links[0][to_v][x].weight; 
                int total_c = to_c+ e_c;
                int total_w = to_w + e_w;
                int e_d_w = dominate_paths_minimum_weight + e_w;
                if(total_c <0 || total_w <0){
                    cout<<"error"<<endl;
                }
                pareto_pair<int, int> sp(total_c, total_w);
                traversal_tuple<int> tri_tmp;
                tri_tmp.v = e_v;
                tri_tmp.cost = sp.first;
                tri_tmp.weight = sp.second;
                tri_tmp.dominate_paths_minimum_weight = e_d_w;
                pq.push(tri_tmp);
            }
        }

        vector<int> sorted_v;
        for(int counter=0; counter<mark[0].occur.size(); counter++){
            int appear_v = mark[0].occur[counter];
            if(mark[0].exist(appear_v) && !pruned[appear_v] && appear_v!=i){
                sorted_v.push_back(appear_v);
            }
        }
        std::sort(sorted_v.begin(), sorted_v.end());
        //fprintf(contrat_file, "%d:", boundary_id[i]);
        fprintf(contrat_file, "%d:", vid2bid[i]);
        for( int counter=0; counter<sorted_v.size(); counter++){
            int vert = sorted_v[counter];
            //vector<pareto_pair<int,int> > &p_list = mark[0][vert];
            vector<pareto_pair<int,int> > p_list;
            for(int x=0; x<mark[0][vert].size(); x++){
                pareto_pair<int, int> pp = mark[0][vert][x];
                int p_size = p_list.size();
                if(p_size>0 && p_list[p_size-1].first <= pp.first &&
                   p_list[p_size-1].second <= pp.dominate_paths_minimum_weight*APPROX){
                    assert(p_list[p_size-1].dominate_paths_minimum_weight > pp.dominate_paths_minimum_weight);
                    p_list[p_size-1].dominate_paths_minimum_weight = pp.dominate_paths_minimum_weight;
                    assert(pp.dominate_paths_minimum_weight*APPROX >= p_list[p_size-1].second);
                    continue;
                }else{
                    p_list.push_back(pp);
                }
            }
            for(int x=0; x<p_list.size(); x++){
                pareto_pair<int, int> pp = p_list[x];
                //int bid_vert = boundary_id[vert];
                link_list_edge forward_tri(vid2bid[vert], pp.second, pp.first);
                forward_tri.dominate_paths_minimum_weight = pp.dominate_paths_minimum_weight;
                fprintf(contrat_file, "%d %d %d %d ", forward_tri.v, forward_tri.weight, forward_tri.cost, 
                        forward_tri.dominate_paths_minimum_weight);
                link_list_edge backward_tri(vid2bid[i], pp.second, pp.first);
                backward_tri.dominate_paths_minimum_weight = pp.dominate_paths_minimum_weight;
                links[0][vid2bid[i]].push_back(forward_tri);
                links[1][vid2bid[vert]].push_back(backward_tri);
                BG_M++;
            }
        }
        fprintf(contrat_file, "-1\n");
    }


    void load_directed_graph(string graph_name, string boundary_graph_name,
                             string partition_name, string cut_name){
        clock_t load_graph_start = clock();
        FILE* file = fopen(graph_name.c_str(), "r");
        FILE* cut_file = fopen(cut_name.c_str(), "r");
        ifstream partition_file(partition_name);
        int ret_val = fscanf(file, "%d", &N);
        if(ret_val <0){ cout<<"read input graph file error"<<endl; exit(0);}
        cout<<N<<" nodes"<<endl;
        
        init_for_G();
        int u, v;
        int weight, cost, dominate_paths_minimum_weight;
        int lines =0;
        double tmp_cost=0;
        M=0;
        for(; fscanf(file, "%d:", &u)==1; lines++){
            for(; fscanf(file, "%d", &v)==1;){
                if(v==-1) break;
                int ret_val = fscanf(file, "%d %d", &weight, &cost);
                if(ret_val <0){cout<<"input file error"<<endl; exit(0);}
                //cout<<u<<"->"<<v<<" "<<weight<<" "<<cost<<endl;
                link_list_edge forward_tri(v, weight, cost);  
                link_list_edge backward_tri(u, weight, cost);
                if(cost > max_cost) max_cost = cost;
                original_links[0][u].push_back(forward_tri);
                original_links[1][v].push_back(backward_tri);
                tmp_cost=  M*1.0/(M+1) * tmp_cost +cost*1.0/(M+1);
                M++;    
            }
        }
        avg_cost= (int)tmp_cost;
        cout<<"Avg cost:"<<avg_cost<<endl;
        cout<<"number of edges: "<<M<<endl;
      
        for(int xx=0; xx<N; xx++){
            pruned[xx]=1;
            partition[xx] =-1;
        }


        cout<<"reading partition file"<<endl;
        string buf;
        string delimt("============================");
        int current_part=-1;
        while(getline(partition_file, buf)){
            buf.erase(buf.find_last_not_of("\n\r\t")+1); 
            if(buf.compare(delimt) ==0 ){
                getline(partition_file, buf);
                buf.erase(buf.find_last_not_of(": \n\r\t")+1); 
                buf = buf.substr(10);
                int partition_num = atoi(buf.c_str());
                //cout<<partition_num<<endl;
                current_part = partition_num;
            }else{
                stringstream ss(buf);
                int node_id;
                double x_axis, y_axis;
                ss>>node_id>>x_axis>>y_axis;
                node_id--;
                partition[node_id] = current_part;
            }
        }

        for(int x=0; x<N; x++){
            if(partition[x]<0) assert(0);
        }

        cout<<"finished reading partition file"<<endl;
        cout<<"number of partitions:"<<(current_part+1)<<endl;

        partition_boundary.re_allocate(current_part+1);
        partition_boundary.m_num = current_part+1;

        iVector<pair<int,int> > cut_vector;
        int cut_start_v, cut_end_v;
        for(lines =0; fscanf(cut_file, "%d %d", &cut_start_v, &cut_end_v)==2; lines++){
            
            pair<int, int> pr;
            if(cut_start_v < cut_end_v){
                pr.first = cut_start_v;
                pr.second = cut_end_v;
            }else{
                pr.first = cut_end_v;
                pr.second = cut_start_v;
            }
            cut_vector.push_back(pr);
        }

        cut_vector.unique();

        for(int x=0; x<cut_vector.size(); x++){
            cut_start_v = cut_vector[x].first;
            cut_end_v = cut_vector[x].second;
            pruned[cut_start_v]=0;
            //pruned[cut_end_v]=0;
            partition_boundary[partition[cut_start_v]].push_back(cut_start_v);
            partition_boundary[partition[cut_end_v]].push_back(cut_start_v);
        }


        for(int x=0; x<partition_boundary.m_num; x++){
            //cout<<"The "<<x<<"'th partition: "<<partition_boundary[x].size()<<endl;
            partition_boundary[x].unique();
        }

        vid2bid.initialize(N);
        int id=0;
        for(int i=0; i<N; i++){
            if(!pruned[i]){
                vid2bid.insert(i, id);
                id++;
            }
        }
        BG_N = id;

        bid2vid.initialize(BG_N);

        for(int i=0; i<vid2bid.occur.size(); i++){
            int x = vid2bid.occur[i];
            if(vid2bid.exist(x)){
                int bid = vid2bid[x];
                bid2vid.insert(bid, x);
            }
        }
        
        
        init_for_BG();

        int pruned_count=0;
        priority_queue<traversal_tuple<int> > pq;
        int pruned_shortcut_count=0;
        FILE* boundary_file = fopen(boundary_graph_name.c_str(), "r");
        int small_n=0;
        ret_val = fscanf(boundary_file, "%d", &small_n);
        if(ret_val <0){ cout<<"read input graph file error"<<endl; exit(0);}



        for(; fscanf(boundary_file, "%d:", &u)==1; ){
            for(; fscanf(boundary_file, "%d", &v)==1;){
                if(v==-1) break;
                //v= boundary_id[v];
                int ret_val = fscanf(boundary_file, "%d %d %d", &weight, &cost,
                                     &dominate_paths_minimum_weight);
                if(ret_val <0){cout<<"input file error"<<endl; exit(0);}
                //cout<<u<<"->"<<v<<" "<<weight<<" "<<cost<<endl;
                link_list_edge forward_tri(v, weight, cost); 
                forward_tri.dominate_paths_minimum_weight = 
                    dominate_paths_minimum_weight;
                link_list_edge backward_tri(u, weight, cost);
                backward_tri.dominate_paths_minimum_weight = dominate_paths_minimum_weight;
                BG_M++;
                if(cost > max_cost) max_cost = cost;
                links[0][u].push_back(forward_tri);
                links[1][v].push_back(backward_tri);
            }
        }
        
        fclose(boundary_file);
        fclose(file);
        fclose(cut_file);
        clock_t load_graph_stop = clock();
        cout<<"overlay graph: N="<<BG_N<<",  M="<<BG_M<<endl;
        printf( "Load Graph Time: %0.2lf\n",(double)(load_graph_stop-load_graph_start)/CLOCKS_PER_SEC );
    }

    void compute_order_by_degree(){
        cout<<"ordering vertices by degree"<<endl;
        vector<pair<int, long> > sort_by_deg;
        for(int i=0; i<BG_N; i++){
            sort_by_deg.push_back(make_pair(i, 
                (links[0][i].size()+1)*(links[1][i].size()+1)));
        }
        sort(sort_by_deg.begin(), sort_by_deg.end(), pair_comp_func_long_bigfirst);
        for(int i=0; i<BG_N; i++){
            level[i] = sort_by_deg[i].first;
        }
        cout<<"ordering finished"<<endl;
    }

    iHeap<double> lambda_pq;
    iMap<double> dist;
    iMap<int> father;
    iMap<int> coverage;
    map<int,int> tree_root;
    void compute_order_by_coverage_heuristic(){
        srand(time(NULL));
        cout<<"ordering vertices by coverage"<<endl;
        cout<<"constructing trees..."<<endl; 
        int available_samples = BG_M*8;
        cout<<"Sample size: "<<available_samples<<endl;

        tree_father.re_allocate(BG_N);
        tree_father.m_num = BG_N;
        tree_children.re_allocate(BG_N);
        tree_children.m_num = BG_N;
        tree_coverage.re_allocate(BG_N);
        tree_coverage.m_num = BG_N;

        int* perm = new int[BG_N + 1];
        lambda_pq.initialize(BG_N);
        dist.initialize(BG_N);
        father.initialize(BG_N);
        coverage.initialize(BG_N);
        for(int i=0;i<BG_N;i++)
        {
            perm[i]=i;
        }
        
        for(int i=BG_N-1; i ; i--){
            int r_id = rand()%(i+1);
            int temp = perm[i];
            perm[i] = perm[r_id];
            perm[r_id] = temp;
        }

        clock_t start_time, end_time;
        start_time = clock();
        //constructing trees
        int tree_id;
        for( int i=0; i<BG_N; i++){
            if(available_samples<1){
                break;
            }
            int v = perm[i];
            tree_id = i;
            tree_root.insert(make_pair(i, v));
            double lambda = ((double) rand() / (RAND_MAX));
            lambda_Dijkstra(v, lambda, tree_id, available_samples);
        }
        cout<<"Sampled "<<tree_id<<" trees"<<endl;

        cout<<"Calculating coverage"<<endl;
        coverage.clean();
        bidij_1d[0].clean();
        for(int i=0; i<BG_N; i++){
            int temp_coverage =0;
            for(auto element: tree_coverage[i]){
                temp_coverage+=element.second;
            }
            coverage.insert(i, temp_coverage);
            bidij_1d[0].insert(i, -temp_coverage);
            
            //assertion passed
            //assert(coverage[i]>=0);
        }

        //ordering vertices;
        cout<<"Ordering vertices"<<endl;
        for(int i=0; i<BG_N; i++){
            //cout<<"ordered "<<i<<" vertices"<<endl;
            int node_to_order = bidij_1d[0].head();
            
            if(coverage[node_to_order]>0){
                bidij_1d[0].insert(node_to_order, -(coverage[node_to_order]) );
            }else{
                bidij_1d[0].insert(node_to_order, -1 *(links[0][node_to_order].size()+1)*(links[1][node_to_order].size()+1));
            }
            
            
            for(;bidij_1d[0].head()!= node_to_order;){
                //cout<<"lazy updating"<<endl;
                node_to_order = bidij_1d[0].head();
                if(coverage[node_to_order]>0){
                    bidij_1d[0].insert(node_to_order, -(coverage[node_to_order]) );
                }else{
                    bidij_1d[0].insert(node_to_order, -1 *(links[0][node_to_order].size()+1)*(links[1][node_to_order].size()+1));
                }
            }


            v2level[node_to_order] = i;
            level[i] = node_to_order;
            bidij_1d[0].pop();
            //cout<<i<<" "<<node_to_order<<" "<<coverage[node_to_order]<<endl; 
            
            //updating father and children's coverage information
            
            //int elimnated_coverage = coverage[node_to_order];

            if(!coverage.exist(node_to_order) ){
               cout<<node_to_order<<endl; 
               assert(coverage.exist(node_to_order));
            }

            /*if(tree_father[node_to_order].size() >0 && tree_children[node_to_order].size() >0){
                cout<<"updating "<<i<<"'th father and children"<<endl;
            }else{
                cout<<"sampled "<<tree_id<<" trees"<<endl;
                int x;
                cin>>x;
                //assert(0);
            }*/


            for(auto v_father : tree_father[node_to_order]){
                //no need to delete edges
                int check_tree_id = v_father.first;
                int elimnated_coverage = tree_coverage[node_to_order][check_tree_id];
                int check_father_id = v_father.second;
                //if(elimnated_coverage< 0) assert(0);
                coverage[check_father_id] -= elimnated_coverage;
                tree_coverage[check_father_id][check_tree_id]-=elimnated_coverage;
                vector<int>& children_v = tree_children[check_father_id][check_tree_id];

                vector<int>::iterator it;
                //assert( (it = find(children_v.begin(), children_v.end(), node_to_order))!= children_v.end() );
                //assert( (it = find(children_v.begin(), children_v.end(), node_to_order))!= children_v.end() );
                if((it = find(children_v.begin(), children_v.end(), node_to_order))!= children_v.end() ){
                    children_v.erase(it);
                }
                //cout<<"tree id: "<<check_tree_id<<" node id: "<<check_father_id<<" root id: "<< tree_root[check_tree_id]<<endl;
                while(tree_father[check_father_id].find(check_tree_id) != tree_father[check_father_id].end()){
                    //cout<<"tracing tree"<<endl; 
                    //cout<<check_tree_id<<" "<<tree_father[check_father_id][check_tree_id]<<endl;
                    check_father_id = tree_father[check_father_id][check_tree_id];
                    /*if(coverage[check_father_id] -  elimnated_coverage <0){
                        cout<<"strange"<<endl;
                        cout<<check_father_id<<" "<<coverage[check_father_id]<<" "<<elimnated_coverage<<endl;
                    }*/
                    tree_coverage[check_father_id][check_tree_id]-=elimnated_coverage;
                    coverage[check_father_id]-=elimnated_coverage;
                    
                    
                    /*if(coverage[check_father_id] <0){ 
                        cout<<check_father_id<<endl;
                        assert(0);
                    }*/
                }
                //cout<<check_father_id<<endl;
               //exit(0);
            }
            tree_father[node_to_order].clear();
            
            for(auto v_children : tree_children[node_to_order]){
                int check_tree_id = v_children.first;
                vector<int>& v_children_list = v_children.second;
                //cout<<"-----listing children------"<<endl;
                vector<int> children_list;
                children_list.clear();
                 for(int k=0; k<v_children_list.size(); k++){
                    int ith_children = v_children_list[k];
                    //cout<<ith_children<<" ";
                    //int self_eliminate_coverage = tree_coverage[ith_children][tree_id];
                    //cout<<self_eliminate_coverage<<endl;
                    children_list.push_back(ith_children);
                 }

                //cout<<"-----traversing tree------"<<endl;
                 for(int k=0; k<children_list.size(); k++){
                    int vx = children_list[k];
                   
                    //assertion passed
                    //assert( tree_father[vx].find(check_tree_id) != tree_father[vx].end());

                    //cout<<vx<<endl;
                    vector<int> & c_list = tree_children[vx][check_tree_id];
                    for(int j=0; j<c_list.size(); j++){
                        children_list.push_back(c_list[j]);
                    }
                    //assert(tree_coverage[vx][check_tree_id] >=0);
                    coverage[vx]-= tree_coverage[vx][check_tree_id];
                    tree_coverage[vx].erase(check_tree_id);
                    tree_children[vx].erase(check_tree_id);
                    tree_father[vx].erase(check_tree_id);
                 }
            }
            
            tree_children[node_to_order].clear();
            //exit(0);
        }

        end_time = clock();
        ordering_time = (double)(end_time-start_time)/CLOCKS_PER_SEC ;
        cout<<"Ordering time: "<<ordering_time<<endl;
        //exit(0);
    }


    void lambda_Dijkstra(int v, double lambda, int tree_id, int& available_samples){
        lambda_pq.clean();
        dist.clean();
        father.clean();
        coverage.clean();
        lambda_pq.insert(v, 0);
        dist.insert(v, 0); 
        father.insert(v, v);
        while(!lambda_pq.empty()){
            int forward_v = lambda_pq.head();
            double lambda_dist = dist[forward_v];
            int father_of_forward_v =  father[forward_v];
            lambda_pq.pop();
            if(forward_v!= v){
                //assert(father_of_forward_v>=0);
                //assert(forward_v!= father_of_forward_v);
                tree_father[forward_v].insert(make_pair(tree_id, father_of_forward_v));
                if(tree_children[father_of_forward_v].find(tree_id) != tree_children[father_of_forward_v].end()){
                    tree_children[father_of_forward_v][tree_id].push_back(forward_v);
                }else{
                   vector<int> children_of_v;
                   children_of_v.push_back(forward_v);
                   tree_children[father_of_forward_v].insert(make_pair(tree_id, children_of_v) );
                }
                available_samples--;
            }
            for(int i=0; i<links[0][forward_v].size(); i++){
                //cout<<"x"<<endl;
                int expand_v = links[0][forward_v][i].v;
                int expand_cost = links[0][forward_v][i].cost;
                int expand_weight = links[0][forward_v][i].weight;
                double expand_lambda_dist = (1-lambda)*expand_cost + lambda*expand_weight;
                double total_lambda_dist = lambda_dist + expand_lambda_dist;
                if(dist.exist(expand_v)){
                    if(total_lambda_dist < dist[expand_v]){
                        dist.insert(expand_v, total_lambda_dist);
                        lambda_pq.insert(expand_v, total_lambda_dist);
                        father.insert(expand_v, forward_v);
                    }
                }else{
                    dist.insert(expand_v, total_lambda_dist);
                    lambda_pq.insert(expand_v, total_lambda_dist);
                    father.insert(expand_v, forward_v);
                }
            }
        }

        vector<pair<double, pair<int, int> > > tree_vec;
        for(int k=0; k<dist.occur.size(); k++){
            int vx = dist.occur[k];
            if(vx!= v && dist.exist(vx) ){
                double v_dist = dist[vx];
                int v_father = father[vx];
                //assert(v_father>=0);
                tree_vec.push_back(make_pair(v_dist, make_pair(vx, v_father) ) ); 
            }
        }
        sort( tree_vec.begin(), tree_vec.end(), lambda_tree_sorter);

        //cout<<tree_vec.size()<<endl;
        for(int i=0; i< tree_vec.size(); i++){
            //cout<<std::setprecision(15) <<tree_vec[i].first<<endl;
            int vx = tree_vec[i].second.first;
            int vx_father = tree_vec[i].second.second;
            if(coverage.exist(vx)){
                coverage[vx]++;    
            }else{
                coverage.insert(vx, 1);
            }
            
            if(vx_father!= vx){
                if(coverage.exist(vx_father)){
                    coverage[vx_father] +=coverage[vx];
                }else{
                    coverage.insert(vx_father, coverage[vx]);
                }
            }
            
            //assertion passed
            //assert(coverage[vx] <=BG_N && coverage[vx_father] <= BG_N);
        }

        for(int i=0; i<coverage.occur.size(); i++){
            int vx= coverage.occur[i];
            if(coverage.exist(vx)){
                //cout<<coverage[vx]<<endl;
                assert(coverage[vx]<= BG_N);
                tree_coverage[vx].insert(make_pair(tree_id, coverage[vx]) );
            }
        }
        //cout<<coverage[ (int) rand()%BG_N]<<endl;
        //exit(0);
    }

    void load_query( string query_file_name){
        FILE* query_file = fopen(query_file_name.c_str(), "r");
        int u, v, cost_lower, cost_upper, cost_limit;
        int cnt=0;
        double tmpcost;
        for(; fscanf(query_file, "%d %d %lf\n", &u, &v, &tmpcost)==3;){
            cost_limit = tmpcost;
            u--;
            v--;
            query_set[cnt].u = u;
            query_set[cnt].v = v;
            //query_set[cnt].cost_limit = cost_limit+ 4*cost_upper;
            query_set[cnt].cost_limit = cost_limit;
            cost_lower = cost_limit;
            cost_upper = cost_limit;
            cnt++;
        }
        fclose(query_file);
    }
    
    int find_out_pos(int weight, vector<triple<int> > &label_list){
        
        int l, r;
        for(l=0,r=label_list.size(); l<r;){
            int m = l+ (r-l)/2;
            if(label_list[m].weight > weight ) l=m+1;
            else r = m;
        }
        return l;
    }

    int find_in_pos(int cost,  vector<triple<int> > & label_list){
        if(cost<0) return label_list.size();
        int l, r;
        for(l=0,r=label_list.size(); l<r;){
            int m = l+ (r-l)/2;
            if(label_list[m].cost > cost ) l=m+1;
            else r = m;
        }
        return l;
    }

    bool opt_is_dominated_by_existing_labelling(int u, int v, int cost, int weight, int dir){
        
        int start_v = dir==0?v2level[u]:v2level[v];
        int end_v = dir==0?v2level[v]:v2level[u];
         label_structure<int>& out_label = opt_labels[0][start_v];
         label_structure<int> & in_label = opt_labels[1][end_v];

         for(int i=0, j=0; i<out_label.label_nodes.size() && j <in_label.label_nodes.size(); ){
            label_cache<int> & lc1 = out_label.label_nodes[i];
            label_cache<int> & lc2 = in_label.label_nodes[j];
            if(lc1.v<lc2.v) i++;
            else if(lc1.v > lc2.v) j++;
            else if( lc1.min_cost + lc2.min_cost > cost) {
                i++;
                j++;
            }else{
               for( int k=lc1.start_pos, m = lc2.start_pos; k<=lc1.end_pos && m <= lc2.end_pos; ){
                   int out_cost = out_label.label_content[k].cost;
                   int out_weight = out_label.label_content[k].weight;

                   int in_cost = in_label.label_content[m].cost;
                   int in_weight = in_label.label_content[m].weight;

                   int total_cost = out_cost + in_cost;
                   int total_weight = out_weight + in_weight;

                   if(total_cost>cost){
                       m++;
                   }else{
                       if(total_weight<= weight){ 
                           return true;
                       }else{
                           k++;
                       }
                   }
               }
               i++;
               j++;
            }
         }
        return false;
    }
    bool is_dominated_by_existing_labelling(int u, int v, int cost, int weight, int dir){
        
        int start_v = dir==0?v2level[u]:v2level[v];
        int end_v = dir==0?v2level[v]:v2level[u];
        for(unsigned int i=0, j=0; i<labels[0][start_v].size() && j<labels[1][end_v].size();){
            int x = label_nodes[0][start_v][i];
            int y = label_nodes[1][end_v][j];
            if( x < y ) i++;
            else if( x > y) j++;
            else{
                vector<triple<int> > & out_list = labels[0][start_v][i];
                vector<triple<int> > & in_list = labels[1][end_v][j];
                if(out_list.size() == 0 || in_list.size() == 0) continue;

                for(unsigned int k=0, m=0/*find_in_pos(cost - out_list[0].cost, in_list)*/; k<out_list.size() && m<in_list.size();){
                    int out_cost = out_list[k].cost;
                    int out_weight = out_list[k].weight;
           

                    //int dominated_out_cost = out_list[k].dominated_cost;
                    //int dominated_out_weight = out_list[k].dominated_weight;

                    int in_cost = in_list[m].cost;
                    int in_weight = in_list[m].weight;
                    int dominated_in_cost = in_list[m].dominated_cost;
                    int dominated_in_weight = in_list[m].dominated_weight;

                    int total_cost = out_cost + in_cost;
                    int total_weight = out_weight + in_weight;
                    //int dominated_total_cost = dominated_out_cost + dominated_in_cost;
                    //int dominated_total_weight = dominated_out_weight + dominated_in_weight;
                    if(total_cost>cost){
                        m++;
                    }else{
                        if(total_weight<= weight){ 
                            return true;
                        }else{
                            k++;
                        }
                    }
                }
                i++;
                j++;
            }
        }

        return false;
    }

    void inline reserve_children_paths_pre_prunning(int v, int t_node, int pl, int wei, int cst, int dir, triple<int>& t,  pareto_pair<int, int>& pareto_pp){
        int pivot_node = level[pl];
        int t_node_level = v2level[t_node];
        path_hash_keys phk_left_child(v, pivot_node, pareto_pp.pivot_first, pareto_pp.pivot_second, dir);
        pre_paths[dir][pivot_node].insert(pair<path_hash_keys, int> (phk_left_child,1));
        if(dir==0){
            t.left_child_first_pointer = pareto_pp.pivot_first;
            t.left_child_second_pointer = pareto_pp.pivot_second;

        }else{

            t.right_child_first_pointer = pareto_pp.pivot_first;
            t.right_child_second_pointer = pareto_pp.pivot_second;
        }
    }

    void node_traverse(int v, int exam_vertex, int s2examv_cost, int s2examv_weight,
                                 int dominate_paths_minimum_weight, int s2examv_pivot_level, int s2examv_pivot_cost, 
                                 int s2examv_pivot_weight,  int dir, int i){
        for(int j=0; j<links[dir][exam_vertex].size(); j++){
            int expand_v = links[dir][exam_vertex][j].v; 
            if(v2level[expand_v]<=i){
                continue;
            }

            int weight = links[dir][exam_vertex][j].weight; 
            int cost = links[dir][exam_vertex][j].cost;
            int check_cost = s2examv_cost + cost; 
            int check_weight = s2examv_weight + weight;
            int extd_domininated_paths_min_weight = links[dir][exam_vertex][j].dominate_paths_minimum_weight; 
            pareto_pair<int, int> skyline_pair(check_cost, check_weight);
            skyline_pair.dominate_paths_minimum_weight = dominate_paths_minimum_weight + extd_domininated_paths_min_weight;
            if(exam_vertex!= v){
                skyline_pair.pivot_level = v2level[exam_vertex];
                skyline_pair.pivot_first = s2examv_cost;
                skyline_pair.pivot_second = s2examv_weight;
            }else{
                skyline_pair.pivot_level = -2;
            }
            
            if(mark[dir].exist(expand_v)){
                vector<pareto_pair<int, int> > &expand_v_list =  mark[dir][expand_v];
                int ev_list_size = expand_v_list.size();
                if(expand_v_list[ev_list_size-1].first<= skyline_pair.first 
                   && expand_v_list[ev_list_size-1].dominate_paths_minimum_weight <= skyline_pair.dominate_paths_minimum_weight){ 
                    continue;
                }
                
                traversal_tuple<int> tri_tmp;
                tri_tmp.v = expand_v;
                tri_tmp.cost = skyline_pair.first;
                tri_tmp.weight = skyline_pair.second;
                tri_tmp.pivot_level = skyline_pair.pivot_level;
                tri_tmp.dominate_paths_minimum_weight = skyline_pair.dominate_paths_minimum_weight;
                tri_tmp.to_pivot_cost = skyline_pair.pivot_first;
                tri_tmp.to_pivot_weight = skyline_pair.pivot_second;
                if(tra_estimate.exist(expand_v) ){
                    if(tra_estimate[expand_v].cost<= tri_tmp.cost && tra_estimate[expand_v].dominate_paths_minimum_weight <= tri_tmp.dominate_paths_minimum_weight){
                        continue;
                    }else{
                        if(tra_estimate[expand_v].cost < tri_tmp.cost)
                            tra_estimate.insert(expand_v, tri_tmp);
                    }
                }else{
                    tra_estimate.insert(expand_v, tri_tmp);
                }

                individual_pq[expand_v].push(tri_tmp);
                //tra_estimate.insert(expand_v, individual_pq[expand_v].top());
                tra_pq.insert(expand_v, individual_pq[expand_v].top());

                //pq.push(tri_tmp);
            }else{
                traversal_tuple<int> tri_tmp;
                tri_tmp.v = expand_v;
                tri_tmp.cost = skyline_pair.first;
                tri_tmp.weight = skyline_pair.second;
                tri_tmp.dominate_paths_minimum_weight = skyline_pair.dominate_paths_minimum_weight;
                tri_tmp.pivot_level = skyline_pair.pivot_level;
                tri_tmp.to_pivot_cost = skyline_pair.pivot_first;
                tri_tmp.to_pivot_weight = skyline_pair.pivot_second;
                
                if(tra_estimate.exist(expand_v) ){
                    if(tra_estimate[expand_v].cost<= tri_tmp.cost && tra_estimate[expand_v].dominate_paths_minimum_weight <= tri_tmp.dominate_paths_minimum_weight){
                        continue;
                    }else{
                        if(tra_estimate[expand_v].cost < tri_tmp.cost)
                            tra_estimate.insert(expand_v, tri_tmp);
                    }
                }else{
                    tra_estimate.insert(expand_v, tri_tmp);
                }
                
                
                individual_pq[expand_v].push(tri_tmp);
                //tra_estimate.insert(expand_v, individual_pq[expand_v].top());
                tra_pq.insert(expand_v, individual_pq[expand_v].top());
                
                //pq.push(tri_tmp);
            }
        }
    }

    void overlay_graph_single_dimension_by_cost(int u, int dir){
       pq_mark.clean();
       pq_mark_max_weight.clean();
       pq_mark.insert(u, 0);
       pq_mark_max_weight.insert(u, 0);
       bidij_1d[dir].clean();
       bidij_1d[dir].insert(u, 0);
       while(!bidij_1d[dir].empty()){
            int forward_v  = bidij_1d[dir].head();
            int cost = pq_mark[forward_v];
            int weight = pq_mark_max_weight[forward_v];
            bidij_1d[dir].pop();
            for(int i=0; i<links[dir][forward_v].size(); i++){
                int expand_v = links[dir][forward_v][i].v;
                if(v2level[expand_v] <= v2level[u] ) continue;
                int expand_cost = links[dir][forward_v][i].cost;
                int expand_weight = links[dir][forward_v][i].weight;
                int total_cost = cost + expand_cost;
                int total_weight = weight + expand_weight;
                if(pq_mark.exist(expand_v)){
                    if(total_cost < pq_mark[expand_v]){
                        pq_mark.insert( expand_v, total_cost);
                        pq_mark_max_weight.insert( expand_v, total_weight);
                        bidij_1d[dir].insert(expand_v, total_cost);
                    }else if (total_cost == pq_mark[expand_v]){
                        if(pq_mark_max_weight[expand_v] > total_weight){
                            pq_mark_max_weight[expand_v] = total_weight;
                        }
                    }
                }else{
                    pq_mark.insert( expand_v, total_cost );
                    pq_mark_max_weight.insert( expand_v, total_weight);
                    bidij_1d[dir].insert(expand_v, total_cost);
                }
            }
       }
    }

    void overlay_graph_single_dimension_by_weight(int u, int dir){
        pq_mark_min_weight.clean();
        pq_mark_min_weight.insert(u, 0);
        bidij_1d[dir].clean();
        bidij_1d[dir].insert(u, 0);
        while(!bidij_1d[dir].empty()){
            int forward_v  = bidij_1d[dir].head();
            int weight = pq_mark_min_weight[forward_v];
            bidij_1d[dir].pop();
            for(int i=0; i<links[dir][forward_v].size(); i++){
                int expand_v = links[dir][forward_v][i].v;
                if(v2level[expand_v] <= v2level[u] ) continue;
                int expand_weight = links[dir][forward_v][i].dominate_paths_minimum_weight;
                int total_weight = weight + expand_weight;
                if(pq_mark_min_weight.exist(expand_v)){
                    if(total_weight < pq_mark_min_weight[expand_v]){
                        pq_mark_min_weight.insert( expand_v, total_weight);
                        bidij_1d[dir].insert(expand_v, total_weight);
                    }
                }else{
                    pq_mark_min_weight.insert( expand_v, total_weight );
                    bidij_1d[dir].insert(expand_v, total_weight);
                }
            }
        }
    }

    void pruned_Dijkstra(int* level){
        for(int i=0; i<BG_N; i++){
            v2level[level[i] ] = i;
        }
        int buf1_help=0;
        int buf2_help=0;
        int no_help=0;
        clock_t start_index_construction = clock();
       
        for(int i=0; i<BG_N; i++){
            int v = level[i];
            //if(pruned[v]) continue;
            /*if(i%100==0){
                cout<<"traversing "<<i<<endl;
            }*/
            for(int dir=0; dir<2; dir++){
                //cout<<"traversing "<<i <<"  "<<dir<<endl;
                overlay_graph_single_dimension_by_cost(v, dir);
                overlay_graph_single_dimension_by_weight(v, dir);
              
                start_prunning.clean();
                for(int x=0; x<pq_mark_min_weight.occur.size(); x++){
                    int xv = pq_mark_min_weight.occur[x];
                    if(pq_mark_min_weight.exist(xv) &&
                       pq_mark_min_weight[xv] >0){
                        int min_w = pq_mark_min_weight[xv];
                        int max_w = pq_mark_max_weight[xv];
                        double val = max_w*1.0/min_w;
                        int number_of_stored_paths = 1+(int)(log(val)/log(APPROX));
                        start_prunning.insert(xv, number_of_stored_paths);
                    }
                }
                

                traversal_tuple<int> tri(v,0,0);
                tri.pivot_level = -2;
                tri.to_pivot_cost=0;
                tri.to_pivot_weight = 0;
                tri.dominate_paths_minimum_weight =0;
                tra_pq.clean();
                tra_estimate.clean();
                tra_pq.insert(v, tri);
                tra_estimate.insert(v, tri);

                vector<pair< pareto_pair<int, int>, pair<int,int> >  > sorted_pareto_paths;
                sorted_pareto_paths.clear();
                while(!tra_pq.empty()){
                    traversal_tuple<int> t = tra_pq.top();
                    tra_pq.pop(); 

                    int exam_vertex = t.v; 
                    int s2examv_weight = t.weight; int s2examv_cost = t.cost;

                    while(!individual_pq[exam_vertex].empty()){
                        traversal_tuple<int> copy_t = individual_pq[exam_vertex].top();

                        if(copy_t.cost >= t.cost && copy_t.dominate_paths_minimum_weight >=
                           t.dominate_paths_minimum_weight) 
                            individual_pq[exam_vertex].pop();
                        else{
                            break;
                        }
                    }

                    if(!individual_pq[exam_vertex].empty()){
                        traversal_tuple<int> copy_t = individual_pq[exam_vertex].top();
                        tra_pq.insert(exam_vertex, copy_t);
                        //tra_estimate.insert(exam_vertex, copy_t);
                    }
                    /*else{
                      tra_estimate.erase(exam_vertex);
                      }*/

                    int s2examv_pivot_level = t.pivot_level;
                    int s2examv_pivot_cost = t.to_pivot_cost;
                    int s2examv_pivot_weight = t.to_pivot_weight;
                    pareto_pair<int, int> checking_pair(s2examv_cost, s2examv_weight);
                    checking_pair.dominate_paths_minimum_weight = t.dominate_paths_minimum_weight;
                    int dominate_paths_minimum_weight = t.dominate_paths_minimum_weight;

                    if(dominate_paths_minimum_weight< 0){
                        assert(0);
                    }

                    if(mark[dir].exist(exam_vertex)){
                        int l_sz = mark[dir][exam_vertex].size();
                        if(mark[dir][exam_vertex][l_sz-1].first <= s2examv_cost && 
                           mark[dir][exam_vertex][l_sz-1].dominate_paths_minimum_weight <=dominate_paths_minimum_weight){
                            continue;
                        }else{
                            if( opt_is_dominated_by_existing_labelling(v, exam_vertex, s2examv_cost, dominate_paths_minimum_weight , dir) ){
                                continue;
                            }

                            checking_pair.pivot_level = t.pivot_level;
                            checking_pair.pivot_first = t.to_pivot_cost;
                            checking_pair.pivot_second = t.to_pivot_weight;
                            if(l_sz+1 > start_prunning[exam_vertex] ){
                                checking_pair.dominate_paths_minimum_weight = max( pq_mark_min_weight[exam_vertex],  int (ceil(s2examv_weight/APPROX) ) );
                            }else{
                                checking_pair.dominate_paths_minimum_weight = dominate_paths_minimum_weight;
                            }
                            mark[dir][exam_vertex].push_back(checking_pair);
                            int pos = mark[dir][exam_vertex].size();
                            if(exam_vertex != v){
                                pair<pareto_pair<int, int>, pair<int, int> > new_path = make_pair(mark[dir][exam_vertex][pos-1], make_pair(exam_vertex,pos-1) );
                                sorted_pareto_paths.push_back(new_path);
                            }
                        }
                    }else{
                        if( opt_is_dominated_by_existing_labelling(v, exam_vertex, s2examv_cost, dominate_paths_minimum_weight , dir) ){
                            continue;
                        }
                        if(start_prunning[exam_vertex] > 1 ){
                            checking_pair.dominate_paths_minimum_weight = max( pq_mark_min_weight[exam_vertex],  int (ceil(s2examv_weight/APPROX) ) );
                        }else{
                            checking_pair.dominate_paths_minimum_weight = dominate_paths_minimum_weight;
                        }
                        checking_pair.pivot_level = t.pivot_level;
                        checking_pair.pivot_first = t.to_pivot_cost;
                        checking_pair.pivot_second = t.to_pivot_weight;
                        checking_pair.direction = -1;
                        vector<pareto_pair<int, int> > p_list;
                        p_list.push_back(checking_pair);
                        mark[dir].insert(exam_vertex, p_list);
                        int pos = mark[dir][exam_vertex].size();
                        if(exam_vertex != v){
                            pair<pareto_pair<int, int>, pair<int, int> > new_path = make_pair(mark[dir][exam_vertex][pos-1], make_pair(exam_vertex,pos-1) );
                            sorted_pareto_paths.push_back(new_path);
                        }
                    }

                    node_traverse(v, exam_vertex, s2examv_cost, s2examv_weight, checking_pair.dominate_paths_minimum_weight,
                                            s2examv_pivot_level, s2examv_pivot_cost, s2examv_pivot_weight,  dir, i);
                }

                //now we add labels
                adding_labels(v, dir, i, sorted_pareto_paths);

                //clean the mark
                if((i+1)%100==0){
                    mark[dir].initialize(N);
                }else{
                    mark[dir].clean();
                }
            }
        }
        
        clock_t end_index_construction = clock();
        //cout<<"total # of labels: "<<total_labels<<", total # of reduced labels: "<<approximate_reduced_labels<<endl;
        //cout<<"compression ratio: "<<(1.0*total_labels/(total_labels+approximate_reduced_labels))<<endl;
        index_construction_time = (double)(end_index_construction-start_index_construction)/CLOCKS_PER_SEC;
        cout<<"Index Construction Time: "<< index_construction_time<<endl;
        cout<<"Preprocessing time: "<<(generate_boundary_time + ordering_time + index_construction_time)<<endl; 
    }

    void adding_labels(int v, int dir, int i, vector<pair< pareto_pair<int, int>, pair<int, int> > >& sorted_pareto_paths){
        triple<int> lb(i, 0,0);
        vector<triple<int> > lb_list;
        lb_list.push_back(lb);
        labels[dir][i].push_back(lb_list);
        label_nodes[dir][i].push_back(i);

        for( unsigned int counter=0; counter<mark[dir].occur.size(); counter++){
            int vx = mark[dir].occur[counter];
            if(mark[dir].exist(vx) && vx != v){
                int pos = mark[dir][vx].size();
                last_pruned[vx].first = mark[dir][vx][pos-1].first;
                last_pruned[vx].second = mark[dir][vx][pos-1].second;
            }
        }
      
        label_mark.clean();
        reverse(sorted_pareto_paths.begin(), sorted_pareto_paths.end());

        int tmp_total=0;
        map<path_hash_keys, int>::iterator it;


        for(int x=0; x< sorted_pareto_paths.size(); x++){
           int t_node = sorted_pareto_paths[x].second.first;
           int position = sorted_pareto_paths[x].second.second;
           pareto_pair<int, int> & tmp_pareto_path = sorted_pareto_paths[x].first;
           int cst = tmp_pareto_path.first;
           int wei = tmp_pareto_path.second;
           int pl = tmp_pareto_path.pivot_level;
           triple<int> t (i, wei, cst);
           t.pivot_level = pl;
           path_hash_keys phk(v, t_node, cst, wei, dir);

           if(position==0 || (pre_paths[dir].exist(t_node) && (it=pre_paths[dir][t_node].find(phk))!=pre_paths[dir][t_node].end()) ){
                //preserved path, should not be pruned 
                if(pl>=0){
                    //we need to add left child and right child to pre_paths and should not prune such paths
                    reserve_children_paths_pre_prunning(v, t_node, pl, wei, cst, dir, t,  tmp_pareto_path);
                }
                

                t.dominated_cost = last_pruned[t_node].first;
                t.dominated_weight = last_pruned[t_node].second;


                if(label_mark.exist(t_node)){
                    label_mark[t_node].push_back(t);
                }else{
                    vector<triple<int> > tmp_vec;
                    tmp_vec.push_back(t);
                    label_mark.insert(t_node, tmp_vec);
                }

                pre_paths[dir][t_node].erase(phk);
                if(pre_paths[dir][t_node].size()==0){
                    pre_paths[dir].erase(t_node);
                }

                if(position>0){
                    last_pruned[t_node].first = mark[dir][t_node][position-1].first;
                    last_pruned[t_node].second = mark[dir][t_node][position-1].dominate_paths_minimum_weight;
                }

            }else{
                if( wei<=last_pruned[t_node].second * SQRT_APPROX && 
                    mark[dir][t_node][position-1].dominate_paths_minimum_weight > last_pruned[t_node].second){
                    //add the path
                    if(pl>=0){
                        reserve_children_paths_pre_prunning(v, t_node, pl, wei, cst, dir, t,  tmp_pareto_path);
                    }
                    
                    t.dominated_cost = last_pruned[t_node].first;
                    t.dominated_weight = last_pruned[t_node].second;

                    if(label_mark.exist(t_node)){
                        label_mark[t_node].push_back(t);
                    }else{
                        vector<triple<int> > tmp_vec;
                        tmp_vec.push_back(t);
                        label_mark.insert(t_node, tmp_vec);
                    }
                        
                    pre_paths[dir][t_node].erase(phk);
                    if(pre_paths[dir][t_node].size()==0){
                        pre_paths[dir].erase(t_node);
                    }
                    
                    //adjust last_pruned
                    last_pruned[t_node].first = mark[dir][t_node][position-1].first;
                    last_pruned[t_node].second = mark[dir][t_node][position-1].second;
                }else{
                    approximate_reduced_labels++;
                    //the path can be pruned, we just omit the path, no need to adjust last_pruned
                }
            }
        }

        for(int counter=0; counter<mark[dir].occur.size(); counter++){
            int x= mark[dir].occur[counter];
                last_pruned[x].first = -1;
                last_pruned[x].second = -1;
        }

        if(1-dir == 0){
            for(int counter = 0; counter < label_mark.occur.size(); counter++){
                int vx = label_mark.occur[counter];
                if(label_mark.exist(vx)){
                    vector<triple<int> > & lb_list = label_mark[vx];

                    if(lb_list.size()>0){
                        label_cache<int> lc;
                        lc.v = i;
                        std::reverse( lb_list.begin(), lb_list.end());
                        lc.min_cost = lb_list[0].cost;
                        labels[0][v2level[vx]].push_back(lb_list);
                        label_nodes[0][v2level[vx]].push_back(i);
                        
                        lc.start_pos = opt_labels[0][v2level[vx]].label_content.size();
                        lc.end_pos = lc.start_pos+lb_list.size()-1;
                        opt_labels[0][v2level[vx]].label_nodes.push_back(lc);
                        opt_labels[0][v2level[vx]].label_content.push_back(lb_list.data(), lb_list.size());
                        tmp_total+= lb_list.size();
                        total_labels+= lb_list.size();
                    }
                }
            }
        }else{
            for(int counter = 0; counter < label_mark.occur.size(); counter++){
                int vx = label_mark.occur[counter];
                if(label_mark.exist(vx)){
                    vector<triple<int> > & lb_list = label_mark[vx];
                    if(lb_list.size()>0){
                        label_cache<int> lc;
                        lc.v = i;
                        lc.start_pos = opt_labels[1][v2level[vx]].label_content.size();
                        lc.end_pos = lc.start_pos+lb_list.size()-1;
                        lc.min_cost = lb_list[lb_list.size()-1].cost;

                        labels[1][v2level[vx]].push_back(lb_list);
                        label_nodes[1][v2level[vx]].push_back(i);
                        opt_labels[1][v2level[vx]].label_nodes.push_back(lc);
                        opt_labels[1][v2level[vx]].label_content.push_back(lb_list.data(), lb_list.size());
                        tmp_total+= lb_list.size();
                        total_labels+= lb_list.size();
                    }
                }
            }
        }
    }


    pair<triple<int> , triple<int> > csp_sketch_on_boundary_graph(int start_v, int end_v, int cost_limit){
        int smallest_weight = INT_MAX;
        triple<int> left_short_cut(-1,-1,-1);
        triple<int> right_short_cut(-1,-1,-1);
        //cout<<"begin to finding the sketch"<<endl;
        for(unsigned int i=0, j=0; i<labels[0][start_v].size() && j<labels[1][end_v].size();){
            int x = label_nodes[0][start_v][i];
            int y = label_nodes[1][end_v][j];
            if( x < y ) i++;
            else if( x > y) j++;
            else{
                vector<triple<int> > & out_list = labels[0][start_v][i];
                vector<triple<int> > & in_list = labels[1][end_v][j];
                for(unsigned int k=0, m=find_in_pos(cost_limit - out_list[0].cost, in_list); k<out_list.size() && m<in_list.size();){
                    int out_cost = out_list[k].cost;
                    int out_weight = out_list[k].weight;
                    int in_cost = in_list[m].cost;
                    int in_weight = in_list[m].weight;
                    int total_cost = out_cost + in_cost;
                    int total_weight = out_weight + in_weight;
                    if(total_cost>cost_limit){
                        m++;
                    }else{
                        if(total_weight<= smallest_weight){ 
                            //assert(out_list[k].v == in_list[m].v); 
                            left_short_cut = out_list[k];
                            right_short_cut = in_list[m];
                            //return true;
                        }
                        k++;
                    }
                }
                i++;
                j++;
            }
        }
        //if(left_short_cut.cost >=0)
        //    cout<<"have found the sketch"<<endl;
        return make_pair(left_short_cut, right_short_cut);
        //return false;
    }
    
    int num_of_pops=0;
    

    vector<edge_list> dijkstra(int u, int v, int cost_limit){
        vector<edge_list> path;
        mark[0].clean();
        traversal_tuple<int> tt(u, 0, 0);
        pq.push(tt);
        //cout<<"wow"<<endl;

        while(!pq.empty()){
            //cout<<pq.size()<<endl;
            traversal_tuple<int> top_t = pq.top();
            pq.pop();
            int to_v = top_t.v;
            int to_c = top_t.cost;
            int to_w = top_t.weight;

            int new_p_pos = 0;

            if(mark[0].exist(to_v)){
                int l_sz = mark[0][to_v].size();
                if(mark[0][to_v][l_sz-1].first <= to_c && mark[0][to_v][l_sz-1].second <=to_w){
                    continue;
                }else{
                    pareto_pair<int, int> pp(to_c, to_w);
                    pp.pivot_level = top_t.pivot_level;
                    pp.pivot_pos = top_t.pivot_pos;
                    pp.pivot_first = top_t.to_pivot_cost;
                    pp.pivot_second = top_t.to_pivot_weight;
                    new_p_pos = mark[0][to_v].size();
                    mark[0][to_v].push_back(pp);
                }
                if( (to_v!= u && !pruned[to_v]) || to_v ==v ) continue;
            }else{
                pareto_pair<int, int> pp(to_c, to_w);
                vector<pareto_pair<int, int> > p_list;
                pp.pivot_level = top_t.pivot_level;
                pp.pivot_pos = top_t.pivot_pos;
                pp.pivot_first = top_t.to_pivot_cost;
                pp.pivot_second = top_t.to_pivot_weight;
                new_p_pos = 0;
                p_list.push_back(pp);
                mark[0].insert(to_v, p_list);
                if( (to_v!= u && !pruned[to_v]) || to_v ==v ) continue;
                //if(!pruned[to_v]) continue;
            }

            for(int x=0; x<original_links[0][to_v].size(); x++){
                int e_v = original_links[0][to_v][x].v;
                if(partition[e_v] != partition[u]) continue;
                int e_c = original_links[0][to_v][x].cost; 
                int e_w = original_links[0][to_v][x].weight; 
                int total_c = to_c+ e_c;
                int total_w = to_w + e_w;
                if(total_c > cost_limit){
                    continue;
                }
                pareto_pair<int, int> sp(total_c, total_w);
                if(mark[0].exist(e_v)){
                    vector<pareto_pair<int, int> > &expand_v_list =  mark[0][e_v];
                    int ev_list_size = expand_v_list.size();
                    if(expand_v_list[ev_list_size-1].first <= sp.first 
                       && expand_v_list[ev_list_size-1].second <=sp.second){
                        continue;
                    }else{
                        traversal_tuple<int> tri_tmp;
                        tri_tmp.v = e_v;
                        tri_tmp.cost = sp.first;
                        tri_tmp.weight = sp.second;
                        tri_tmp.pivot_pos = new_p_pos;
                        tri_tmp.pivot_level = to_v;
                        tri_tmp.to_pivot_cost = e_c;
                        tri_tmp.to_pivot_weight = e_w;
                        pq.push(tri_tmp);
                    }
                }else{
                    traversal_tuple<int> tri_tmp;
                    tri_tmp.v = e_v;
                    tri_tmp.cost = sp.first;
                    tri_tmp.weight = sp.second;
                    tri_tmp.pivot_pos = new_p_pos;
                    tri_tmp.pivot_level = to_v;
                    tri_tmp.to_pivot_cost = e_c;
                    tri_tmp.to_pivot_weight = e_w;
                    pq.push(tri_tmp);
                }
            }
        }

        vector<edge_list> tmp_path;
    
        int l=-1;
        
        for(int i=0; i<mark[0][v].size(); i++){
            if(mark[0][v][i].first <= cost_limit  ) l=i;
            else break;
        }

        if(l<0) return path;

        reachable_queries++;

        int pre_n = mark[0][v][l].pivot_level;
        int pre_pos = mark[0][v][l].pivot_pos;
        //cout<<pre_n<<" "<<pre_pos<<endl;
        edge_list e(mark[0][v][l].pivot_level, v, mark[0][v][l].pivot_first, mark[0][v][l].pivot_second);
        tmp_path.push_back(e);
        while(pre_n != u){
            pareto_pair<int, int> & pp = mark[0][pre_n][pre_pos];
            edge_list tmp(pp.pivot_level, pre_n, pp.pivot_first, pp.pivot_second);
            tmp_path.push_back(tmp);
            pre_n = pp.pivot_level;
            pre_pos = pp.pivot_pos;
        }

        reverse(tmp_path.begin(), tmp_path.end() );
        path.insert(path.end(), tmp_path.begin(), tmp_path.end());
        return path;
    }
   


    void partition_traversal(int node, int dir, int cost_limit){
        start_prunning.clean();
        in_boundary_single_dimension_by_weight(node, dir, cost_limit);
        in_boundary_single_dimension_by_cost(node, dir, cost_limit);
        for(int k=0; k<pq_mark_min_weight.occur.size(); k++){
            int reachable_v = pq_mark_min_weight.occur[k];
            if(pq_mark_min_weight.exist(reachable_v)){
                int min_w = pq_mark_min_weight[reachable_v];
                int max_w = pq_mark_max_weight[reachable_v];
                if(min_w>0 && max_w >0){
                    double val = max_w*1.0/min_w;
                    int number_of_stored_paths = 1+(int)(log(val)/log(APPROX));
                    start_prunning.insert(reachable_v, number_of_stored_paths);
                }
            }
        }
        mark[dir].clean();
        traversal_tuple<int> tt(node, 0, 0);
        tt.dominate_paths_minimum_weight = 0;
        pq.push(tt);
        while(!pq.empty()){
            traversal_tuple<int> top_t = pq.top();
            pq.pop();
            int to_v = top_t.v;
            int to_c = top_t.cost;
            int to_w = top_t.weight;
            int dominate_paths_minimum_weight = top_t.dominate_paths_minimum_weight;
            int new_p_pos = 0;

            pareto_pair<int, int> pp(to_c, to_w);
            if(mark[dir].exist(to_v)){
                int l_sz = mark[dir][to_v].size();
                if(mark[dir][to_v][l_sz-1].first <= to_c && 
                   mark[dir][to_v][l_sz-1].dominate_paths_minimum_weight <=
                   dominate_paths_minimum_weight){
                    continue;
                }else{
                    pp.pivot_level = top_t.pivot_level;
                    pp.pivot_pos = top_t.pivot_pos;
                    pp.pivot_first = top_t.to_pivot_cost;
                    pp.pivot_second = top_t.to_pivot_weight;
                    new_p_pos = mark[dir][to_v].size();
                    if(l_sz+1 > start_prunning[to_v]){
                        pp.dominate_paths_minimum_weight = max( pq_mark_min_weight[to_v],  int (ceil(to_w/BOUNDARY_GRAPH_APPROX) ) );
                    }else{
                        pp.dominate_paths_minimum_weight = dominate_paths_minimum_weight;
                    }
                    mark[dir][to_v].push_back(pp);
                }
                if(!pruned[to_v] && to_v != node) continue;
            }else{
                pp.pivot_level = top_t.pivot_level;
                pp.pivot_pos = top_t.pivot_pos;
                pp.pivot_first = top_t.to_pivot_cost;
                pp.pivot_second = top_t.to_pivot_weight;
                new_p_pos = 0;
                if(start_prunning[to_v] > 1){
                    pp.dominate_paths_minimum_weight = max( pq_mark_min_weight[to_v],  int (ceil(to_w/BOUNDARY_GRAPH_APPROX) ) );
                }else{
                    pp.dominate_paths_minimum_weight = dominate_paths_minimum_weight;
                }
                vector<pareto_pair<int, int> > p_list;
                p_list.push_back(pp);
                mark[dir].insert(to_v, p_list);
                if(!pruned[to_v] && to_v != node) continue;
            }

            dominate_paths_minimum_weight = pp.dominate_paths_minimum_weight;
            
            for(int x=0; x<original_links[dir][to_v].size(); x++){
                int e_v = original_links[dir][to_v][x].v;
                //if(partition[e_v] != partition[node]) continue;
                int e_c = original_links[dir][to_v][x].cost; 
                int e_w = original_links[dir][to_v][x].weight; 
                int total_c = to_c+ e_c;
                int total_w = to_w + e_w;
                int e_d_w = dominate_paths_minimum_weight + e_w;
                if(total_c > cost_limit){
                    continue;
                }
                traversal_tuple<int> tri_tmp;
                tri_tmp.v = e_v;
                tri_tmp.cost = total_c;
                tri_tmp.weight = total_w;
                tri_tmp.pivot_pos = new_p_pos;
                tri_tmp.dominate_paths_minimum_weight = e_d_w;
                tri_tmp.pivot_level = to_v;
                tri_tmp.to_pivot_cost = e_c;
                tri_tmp.to_pivot_weight = e_w;
                pq.push(tri_tmp);
            }
        }
    }

    int unfold_shortcut_num=0;
    int unfold_edge_num=0;

    int reachable_queries=0;

    vector<edge_list> cross_partition_query(int u, int v, int cost_limit){
        vector<edge_list> path;
       
        //step 1.1, forward dijkstra from u to its partition nodes
        
        partition_traversal(u, 0, cost_limit);
        partition_traversal(v, 1, cost_limit);

        //step 1.2, assemble it with boundary nodes 
        pq_mark.clean();
       
        pq_mark_back.clean();

        vector<int> u_to_nodes;
        vector<vector<triple<int> > > u_to_label;
        //cout<<"out lables"<<endl;
        for(int i=0; i<mark[0].occur.size(); i++){
            int vx = mark[0].occur[i];
            //cout<<vx<<endl;
            if(mark[0].exist(vx)&& vid2bid.exist(vx)){
                int level_vx = v2level[vid2bid[vx]];
                vector<int> &label_vx_nodes = label_nodes[0][level_vx];
                for( int j=0; j<label_vx_nodes.size(); j++){
                    if(!pq_mark_back.exist(label_vx_nodes[j])){
                        u_to_nodes.push_back(label_vx_nodes[j]);
                        pq_mark_back.insert(label_vx_nodes[j], 1);
                        cost_min_query[label_vx_nodes[j]] = mark[0][vx][0].first +
                            labels[0][level_vx][j][0].cost;
                        weight_max_query[label_vx_nodes[j]] = mark[0][vx][0].second
                            + labels[0][level_vx][j][0].weight;
                        int sz = mark[0][vx].size();
                        int sz2 = labels[0][level_vx][j].size();
                        weight_min_query[label_vx_nodes[j]] = mark[0][vx][sz-1].second +
                            (int) (labels[0][level_vx][j][sz2-1].weight/APPROX);

                        //assert(weight_max_query[label_vx_nodes[j]] >= weight_min_query[label_vx_nodes[j]]);
                    }else{
                        if(cost_min_query[label_vx_nodes[j]] > mark[0][vx][0].first +
                            labels[0][level_vx][j][0].cost){
                            cost_min_query[label_vx_nodes[j]] = mark[0][vx][0].first +
                                labels[0][level_vx][j][0].cost;
                            weight_max_query[label_vx_nodes[j]] = mark[0][vx][0].second
                                + labels[0][level_vx][j][0].weight;
                        }

                        int sz = mark[0][vx].size();
                        int sz2 = labels[0][level_vx][j].size();
                        if(weight_min_query[label_vx_nodes[j]] > mark[0][vx][sz-1].second +
                            (int)(labels[0][level_vx][j][sz2-1].weight/APPROX) ){
                            weight_min_query[label_vx_nodes[j]] = mark[0][vx][sz-1].second +
                                (int) (labels[0][level_vx][j][sz2-1].weight/APPROX) ;
                        }
                        //assert(weight_max_query[label_vx_nodes[j]] >= weight_min_query[label_vx_nodes[j]]);
                        
                    }
                }
            }
        }
        
        sort(u_to_nodes.begin(), u_to_nodes.end());
        
        

        //cout<<"in lables"<<endl;
        pq_mark_back.clean();
        vector<int> to_v_nodes;
        vector<vector<triple<int> > > to_v_label;
        for(int i=0; i<mark[1].occur.size(); i++){
            int vx = mark[1].occur[i];
            if(mark[1].exist(vx)&& vid2bid.exist(vx) ){
                int level_vx = v2level[vid2bid[vx]];

                vector<int> &label_vx_nodes = label_nodes[1][level_vx];
                for( int j=0; j<label_vx_nodes.size(); j++){
                    if(!pq_mark_back.exist(label_vx_nodes[j])){
                       to_v_nodes.push_back(label_vx_nodes[j]);
                       pq_mark_back.insert(label_vx_nodes[j], 1);
                        int sz = mark[1][vx].size();
                        int sz2 = labels[1][level_vx][j].size();
                        
                        cost_min_query_2[label_vx_nodes[j]] = mark[1][vx][0].first +
                            labels[1][level_vx][j][sz2-1].cost;
                        weight_max_query_2[label_vx_nodes[j]] = mark[1][vx][0].second
                            + labels[1][level_vx][j][sz2-1].weight;

                        weight_min_query_2[label_vx_nodes[j]] = mark[1][vx][sz-1].second +
                            (int)(labels[1][level_vx][j][0].weight/APPROX);

                        //assert(weight_max_query_2[label_vx_nodes[j]] >= weight_min_query_2[label_vx_nodes[j]]);
                    }else{
                        int sz = mark[1][vx].size();
                        int sz2 = labels[1][level_vx][j].size();
                        if(cost_min_query_2[label_vx_nodes[j]] > mark[1][vx][0].first +
                            labels[1][level_vx][j][sz2-1].cost){
                            cost_min_query_2[label_vx_nodes[j]] = mark[1][vx][0].first +
                                labels[1][level_vx][j][sz2-1].cost;
                            weight_max_query_2[label_vx_nodes[j]] = mark[1][vx][0].second
                                + labels[1][level_vx][j][sz2-1].weight;
                        }

                        if(weight_min_query_2[label_vx_nodes[j]] > mark[1][vx][sz-1].second +
                            (int)(labels[1][level_vx][j][0].weight/APPROX)){
                            weight_min_query_2[label_vx_nodes[j]] = mark[1][vx][sz-1].second +
                                (int)(labels[1][level_vx][j][0].weight/APPROX);
                        }
                        //assert(weight_max_query_2[label_vx_nodes[j]] >= weight_min_query_2[label_vx_nodes[j]]);
                    }

                }
            }
        }

        sort(to_v_nodes.begin(), to_v_nodes.end());

        int weight_upper_bound=INT_MAX;
        int num_common_intermediate=0;

        vector<int> init_common_intermediate;
        for(int x=0,y=0; x<u_to_nodes.size() && y<to_v_nodes.size(); ){

            if(u_to_nodes[x]< to_v_nodes[y]) x++;
            else if( u_to_nodes[x] > to_v_nodes[y]) y++;
            else{
                int tmp_cost_lower_bound = cost_min_query[u_to_nodes[x]]
                    + cost_min_query_2[u_to_nodes[x]];
                int tmp_weight_upper_bound = weight_max_query[u_to_nodes[x]]
                    + weight_max_query_2[u_to_nodes[x]];
                if(tmp_cost_lower_bound>cost_limit)
                {
                    x++; y++;
                    continue;
                }
                if(tmp_weight_upper_bound < weight_upper_bound){ 
                    weight_upper_bound = tmp_weight_upper_bound;
                }
                init_common_intermediate.push_back(u_to_nodes[x]);
                
                x++;
                y++;
            } 
        }
       
        //cout<<"before: "<<init_common_intermediate.size()<<endl;

        //cout<<weight_upper_bound<<endl;

        vector<int> common_intermediate;
        for(int x = 0; x<init_common_intermediate.size(); x++){
            int vx = init_common_intermediate[x];
            if(weight_min_query[vx]+weight_min_query_2[vx] > weight_upper_bound){
                /*assert(weight_max_query[vx] + weight_max_query_2[vx] 
                       >= weight_min_query[vx]+ weight_min_query_2[vx]);
                cout<<weight_max_query[vx]<<" "<<weight_max_query_2[vx]<<" "<<weight_upper_bound<<endl;*/
                continue;
            }
            pq_mark.insert(vx, num_common_intermediate);
            vector<triple<int> >  tl;
            u_to_label.push_back(tl);
            to_v_label.push_back(tl);
            common_intermediate.push_back(vx);
            num_common_intermediate++;
        }
        
       
        //cout<<"After: "<<common_intermediate.size()<<endl;
       
        //reachable_queries++;
       
        int u_partition = partition[u];
        for(int i=0; i<partition_boundary[u_partition].size(); i++){
            int vx = partition_boundary[u_partition][i];
            if(mark[0].exist(vx)){
                int level_vx = v2level[vid2bid[vx]];
                vector<pareto_pair<int, int> > &p_list = mark[0][vx];
                vector<vector<triple<int> > > &label_vx = labels[0][level_vx];
                for( int x=0; x< p_list.size(); x++){
                    for(int y=0; y<label_vx.size(); y++){
                        if(label_vx[y].size()>0){
                            int u_to_l = label_vx[y][0].v;
                            if(!pq_mark.exist(u_to_l)) continue;
                            int pos = pq_mark[u_to_l];
                            for (int z=0; z<label_vx[y].size(); z++){
                                triple<int> tr = label_vx[y][z];
                                tr.cost = tr.cost + p_list[x].first;
                                /*if(tr.cost<= cost_limit){
                                }else break;
                                if(x<p_list.size() -1 && p_list[x+1].first + label_vx[y][0].cost <= tr.cost 
                                   && p_list[x+1].second + label_vx[y][0].weight <= tr.weight){
                                    break;
                                }*/
                                tr.weight = tr.weight + p_list[x].second;
                                tr.left_child_first_pointer = vx;
                                tr.left_child_second_pointer = x;
                                tr.pivot_level = level_vx;
                                tr.right_child_first_pointer = y;
                                tr.right_child_second_pointer = z; 
                                u_to_label[pos].push_back(tr);
                            }
                        }
                    }
                }
            }
        }

        for(int i=0; i< u_to_label.size(); i++){
            sort(u_to_label[i].begin(), u_to_label[i].end(), triple_small_first);
            //u_to_label[i].erase(unique(u_to_label[i].begin(), u_to_label[i].end()), u_to_label[i].end()); 
            vector<triple<int> > ll;
            for(int j=0; j<u_to_label[i].size(); j++){
                if(ll.size()>0 && ll[ll.size()-1].cost <= u_to_label[i][j].cost &&
                   ll[ll.size()-1].weight <= u_to_label[i][j].weight) 
                    continue;
                ll.push_back(u_to_label[i][j]);
            }
            u_to_label[i] = ll;
        }

        int v_partition = partition[v];
        for(int i=0; i<partition_boundary[v_partition].size(); i++){
            int vx = partition_boundary[v_partition][i];
            if(mark[1].exist(vx)){
                int level_vx = v2level[vid2bid[vx]];
                vector<pareto_pair<int, int> > &p_list = mark[1][vx];
                vector<vector<triple<int> > > &label_vx = labels[1][level_vx];
                for( int x=0; x< p_list.size(); x++){
                    for(int y=0; y<label_vx.size(); y++){
                        if(label_vx[y].size()>0){
                            int to_v_l = label_vx[y][0].v;
                            if(!pq_mark.exist(to_v_l)) continue;
                            int pos = pq_mark[to_v_l];
                            for (int z=0; z<label_vx[y].size(); z++){
                                triple<int> tr = label_vx[y][z];
                                tr.cost = tr.cost + p_list[x].first;
                                //if(tr.cost> cost_limit) break;
                                //cout<<p_list.size()<<" "<<label_vx[y].size()<<endl;
                                /*if(x<p_list.size() -1 && p_list[x+1].first + label_vx[y][0].cost <= tr.cost 
                                   && p_list[x+1].second + label_vx[y][0].weight <= tr.weight){
                                    cout<<p_list.size()<<" "<<label_vx[y].size()<<endl;
                                    break;
                                }*/
                                tr.right_child_first_pointer = vx;
                                tr.right_child_second_pointer = x;
                                tr.pivot_level = level_vx;
                                tr.left_child_first_pointer = y;
                                tr.left_child_second_pointer = z; 
                                tr.weight = tr.weight + p_list[x].second;
                                to_v_label[pos].push_back(tr);
                            }
                        }
                    }
                }
            }
        }
        

        for(int i=0; i< to_v_label.size(); i++){
            sort(to_v_label[i].begin(), to_v_label[i].end(), triple_small_first);
            vector<triple<int> > ll;
            for(int j=0; j<to_v_label[i].size(); j++){
                if(ll.size()>0 && ll[ll.size()-1].cost <= to_v_label[i][j].cost &&
                   ll[ll.size()-1].weight <=to_v_label[i][j].weight){
                    continue;
                }
                ll.push_back(to_v_label[i][j]);
            }
            reverse(ll.begin(), ll.end());
            to_v_label[i] = ll;
        }

        //step 3 lable intersection and retrieve the constraint shortest path sketch 

        triple<int> left_short_cut(-1,-1,-1);
        triple<int> right_short_cut(-1,-1,-1);

        int smallest_weight = INT_MAX;

        for( int i=0; i< common_intermediate.size();i++ ){
            vector<triple<int> > & out_list = u_to_label[i];
            vector<triple<int> > & in_list = to_v_label[i];
            for(unsigned int k=0, m=0 ; k<out_list.size() && m<in_list.size();){
                //assert(out_list[k].v == in_list[m].v);
                int out_cost = out_list[k].cost;
                int out_weight = out_list[k].weight;
                int in_cost = in_list[m].cost;
                int in_weight = in_list[m].weight;
                int total_cost = out_cost + in_cost;
                int total_weight = out_weight + in_weight;
                if(total_cost>cost_limit){
                    m++;
                }else{
                    if(total_weight<= smallest_weight){ 
                        left_short_cut = out_list[k];
                        right_short_cut = in_list[m];
                    }
                    k++;
                }
            }
        }

        //step 4 unfolding the path
        if(left_short_cut.cost < 0 || right_short_cut.cost <0) return path;

        reachable_queries++;
        //cout<<"cost_limit:"<<cost_limit<<" cost:"<<(left_short_cut.cost+right_short_cut.cost)<<" weight:"<<(left_short_cut.weight+ right_short_cut.weight)<<endl;

        //cout<<u<<" -> "<<left_short_cut.left_child_first_pointer<<" | ";
        pareto_pair<int, int> left_boundary_path = mark[0][left_short_cut.left_child_first_pointer][left_short_cut.left_child_second_pointer];

        int pre_n = left_boundary_path.pivot_level;
        int pre_pos = left_boundary_path.pivot_pos;
        if(pre_n >=0){
            vector<edge_list> tmp_path;
            //cout<<pre_n<<" "<<pre_pos<<endl;
            edge_list e( left_boundary_path.pivot_level, left_short_cut.left_child_first_pointer, left_boundary_path.pivot_first, left_boundary_path.pivot_second);
            tmp_path.push_back(e);
            while(pre_n != u){
                pareto_pair<int, int> & pp = mark[0][pre_n][pre_pos];
                edge_list tmp(pp.pivot_level, pre_n, pp.pivot_first, pp.pivot_second);
                tmp_path.push_back(tmp);
                pre_n = pp.pivot_level;
                pre_pos = pp.pivot_pos;
            }
            reverse(tmp_path.begin(), tmp_path.end() );
            path.insert(path.end(), tmp_path.begin(), tmp_path.end());
        }



        //assert(left_short_cut.left_child_first_pointer == level[left_short_cut.pivot_level]);

        triple<int> left_link = labels[0][left_short_cut.pivot_level][left_short_cut.right_child_first_pointer][left_short_cut.right_child_second_pointer];
        //cout<<level[left_short_cut.pivot_level]<< " -> "<<level[left_link.v]<<" | ";

        triple<int> right_link = labels[1][right_short_cut.pivot_level][right_short_cut.left_child_first_pointer][right_short_cut.left_child_second_pointer];
        //cout<<level[right_link.v]<< " -> "<<level[right_short_cut.pivot_level]<<" | ";

        //cout<<right_short_cut.right_child_first_pointer<< " ->  "<<v<<endl;
        pareto_pair<int, int> right_boundary_path = mark[1][right_short_cut.right_child_first_pointer][right_short_cut.right_child_second_pointer];

        vector<edge_list> shortcut_unfold_list_left;
        vector<edge_list> shortcut_unfold_list_right;


        //cout<<left_short_cut.pivot_level<<" "<<left_link.pivot_level<<" "<<left_link.v<<endl;
        path_unfolding(left_link, left_short_cut.pivot_level, left_link.v, 0, shortcut_unfold_list_left );
      
        //cout<<"unfolding left list boundary edges"<<endl;
        for(int k=0; k<shortcut_unfold_list_left.size(); k++){
            //cout<<v2level[shortcut_unfold_list_left[k].u]<<" "<< v2level[shortcut_unfold_list_left[k].v]<<" "<<
            //    shortcut_unfold_list_left[k].cost<<" "<<shortcut_unfold_list_left[k].weight<<endl;
            if(partition[shortcut_unfold_list_left[k].u] == partition[shortcut_unfold_list_left[k].v]){
                unfold_shortcut_num++;
                shortcut_unfold(bid2vid[shortcut_unfold_list_left[k].u], bid2vid[shortcut_unfold_list_left[k].v], 
                                shortcut_unfold_list_left[k].cost, shortcut_unfold_list_left[k].weight, path);
            }else{
                path.push_back(shortcut_unfold_list_left[k]);
                unfold_edge_num++;
            }
        }

        return path;
        //cout<<"finished unfold left path"<<endl;
        
        //cout<<right_link.v<<" "<<right_link.pivot_level<<" "<<right_short_cut.pivot_level<<endl;
        path_unfolding(right_link, right_link.v, right_short_cut.pivot_level, 1, shortcut_unfold_list_right);
       


        //cout<<"unfolding right list boundary edges"<<endl;
        for(int k=0; k<shortcut_unfold_list_right.size(); k++){
            //cout<<v2level[shortcut_unfold_list_right[k].u]<<" "<< v2level[shortcut_unfold_list_right[k].v]<<" "<<
            //    shortcut_unfold_list_right[k].cost<<" "<<shortcut_unfold_list_right[k].weight<<endl;
            if(partition[shortcut_unfold_list_right[k].u] == partition[shortcut_unfold_list_right[k].v]){
                shortcut_unfold(bid2vid[shortcut_unfold_list_right[k].u], bid2vid[shortcut_unfold_list_right[k].v], 
                                shortcut_unfold_list_right[k].cost, shortcut_unfold_list_right[k].weight, path);
                unfold_shortcut_num++;
            }else{
                path.push_back(shortcut_unfold_list_right[k]);
                unfold_edge_num++;
            }
        }

        pre_n = right_boundary_path.pivot_level;
        pre_pos = right_boundary_path.pivot_pos;
        if(pre_n >=0){
           vector<edge_list> tmp_path;
            //cout<<pre_n<<" "<<pre_pos<<endl;
            edge_list e( right_short_cut.right_child_first_pointer, right_boundary_path.pivot_level, right_boundary_path.pivot_first, right_boundary_path.pivot_second);
            tmp_path.push_back(e);
            while(pre_n != v){
                pareto_pair<int, int> & pp = mark[1][pre_n][pre_pos];
                edge_list tmp(pre_n, pp.pivot_level,  pp.pivot_first, pp.pivot_second);
                path.push_back(tmp);
                pre_n = pp.pivot_level;
                pre_pos = pp.pivot_pos;
            }
        }
        
        return path;
    }


    void shortcut_unfold(int u, int v, int cost_limit, int weight_limit, vector<edge_list> & path){
        mark[0].clean();
        traversal_tuple<int> tt(u, 0, 0);
        pq.push(tt);
        //cout<<"wow"<<endl;

        while(!pq.empty()){
            //cout<<pq.size()<<endl;
            traversal_tuple<int> top_t = pq.top();
            pq.pop();
            int to_v = top_t.v;
            int to_c = top_t.cost;
            int to_w = top_t.weight;

            int new_p_pos = 0;

            if(mark[0].exist(to_v)){
                int l_sz = mark[0][to_v].size();
                if(mark[0][to_v][l_sz-1].first <= to_c && mark[0][to_v][l_sz-1].second <=to_w){
                    continue;
                }else{
                    pareto_pair<int, int> pp(to_c, to_w);
                    pp.pivot_level = top_t.pivot_level;
                    pp.pivot_pos = top_t.pivot_pos;
                    pp.pivot_first = top_t.to_pivot_cost;
                    pp.pivot_second = top_t.to_pivot_weight;
                    new_p_pos = mark[0][to_v].size();
                    mark[0][to_v].push_back(pp);
                }
                if( (to_v!= u && !pruned[to_v]) || to_v ==v ) continue;
            }else{
                pareto_pair<int, int> pp(to_c, to_w);
                vector<pareto_pair<int, int> > p_list;
                pp.pivot_level = top_t.pivot_level;
                pp.pivot_pos = top_t.pivot_pos;
                pp.pivot_first = top_t.to_pivot_cost;
                pp.pivot_second = top_t.to_pivot_weight;
                new_p_pos = 0;
                p_list.push_back(pp);
                mark[0].insert(to_v, p_list);
                if( (to_v!= u && !pruned[to_v]) || to_v ==v ) continue;
                //if(!pruned[to_v]) continue;
            }

            for(int x=0; x<original_links[0][to_v].size(); x++){
                int e_v = original_links[0][to_v][x].v;
                //if(partition[e_v] != partition[u]) continue;
                int e_c = original_links[0][to_v][x].cost; 
                int e_w = original_links[0][to_v][x].weight; 
                int total_c = to_c+ e_c;
                int total_w = to_w + e_w;
                if(total_c > cost_limit || total_w > weight_limit){
                    continue;
                }
                pareto_pair<int, int> sp(total_c, total_w);
                if(mark[0].exist(e_v)){
                    vector<pareto_pair<int, int> > &expand_v_list =  mark[0][e_v];
                    int ev_list_size = expand_v_list.size();
                    if(expand_v_list[ev_list_size-1].first <= sp.first 
                       && expand_v_list[ev_list_size-1].second <=sp.second){
                        continue;
                    }else{
                        traversal_tuple<int> tri_tmp;
                        tri_tmp.v = e_v;
                        tri_tmp.cost = sp.first;
                        tri_tmp.weight = sp.second;
                        tri_tmp.pivot_pos = new_p_pos;
                        tri_tmp.pivot_level = to_v;
                        tri_tmp.to_pivot_cost = e_c;
                        tri_tmp.to_pivot_weight = e_w;
                        pq.push(tri_tmp);
                    }
                }else{
                    traversal_tuple<int> tri_tmp;
                    tri_tmp.v = e_v;
                    tri_tmp.cost = sp.first;
                    tri_tmp.weight = sp.second;
                    tri_tmp.pivot_pos = new_p_pos;
                    tri_tmp.pivot_level = to_v;
                    tri_tmp.to_pivot_cost = e_c;
                    tri_tmp.to_pivot_weight = e_w;
                    pq.push(tri_tmp);
                }
            }
        }

        bool exist_such_path = false;
        vector<edge_list> tmp_path;
    
        int l,r;
        pareto_pair<int, int> sp(cost_limit, weight_limit);
        for(l=0,r=mark[0][v].size(); l<r;){
            int m = l+ (r-l)/2;
            if(mark[0][v][m] < sp ) l=m+1;
            else r = m;
        }

        //assert(l<= mark[0][v].size()-1); 


        if(l>=mark[0][v].size()) return;
        int pre_n = mark[0][v][l].pivot_level;
        int pre_pos = mark[0][v][l].pivot_pos;
        //cout<<pre_n<<" "<<pre_pos<<endl;
        edge_list e(mark[0][v][l].pivot_level, v, mark[0][v][l].pivot_first, mark[0][v][l].pivot_second);
        tmp_path.push_back(e);
        while(pre_n != u){
            pareto_pair<int, int> & pp = mark[0][pre_n][pre_pos];
            edge_list tmp(pp.pivot_level, pre_n, pp.pivot_first, pp.pivot_second);
            tmp_path.push_back(tmp);
            pre_n = pp.pivot_level;
            pre_pos = pp.pivot_pos;
        }

        reverse(tmp_path.begin(), tmp_path.end() );
        path.insert(path.end(), tmp_path.begin(), tmp_path.end());
    }


    void path_unfolding(triple<int> sketch, int start_v, int end_v, int dir, vector<edge_list> & path){
        triple<int> current_sketch = sketch;
        triple<int> next_sketch;
        //cout<<"-------------------------------------"<<endl;
        //cout<<start_v<<" "<<end_v<<" "<<dir<<endl; 

        if(start_v == end_v)
            return;
        
        if(dir==0){
            int current_start =  start_v;
            while(current_sketch.pivot_level >=0 ){
                int node_lv = current_sketch.pivot_level;
                int child_first_pointer = current_sketch.right_child_first_pointer;
                int child_second_pointer = current_sketch.right_child_second_pointer;
                //cout<<dir<<endl;
                //cout<<node_lv<<" "<<child_first_pointer<<" "<<child_second_pointer<<" "<<current_sketch.left_child_first_pointer<<" "<<
                //    current_sketch.left_child_second_pointer<<endl;
                next_sketch = labels[dir][node_lv][child_first_pointer][child_second_pointer];
                edge_list e( level[current_start], level[node_lv], current_sketch.cost - next_sketch.cost, 
                             current_sketch.weight - next_sketch.weight);
                
                path.push_back(e);
                current_start = node_lv;
                current_sketch = next_sketch;
            }
            
            if(current_start != end_v){
                edge_list e ( level[current_start], level[end_v], current_sketch.cost, current_sketch.weight); 
                path.push_back(e);
            }
        }else{
            vector<edge_list>  tmp_path;
            int current_end = end_v;
            while(current_sketch.pivot_level >= 0){
                int node_lv = current_sketch.pivot_level;
                int child_first_pointer = current_sketch.left_child_first_pointer;
                int child_second_pointer = current_sketch.left_child_second_pointer;
                //cout<<dir<<endl;
                //cout<<node_lv<<" "<<child_first_pointer<<" "<<child_second_pointer<<" "<<current_sketch.right_child_first_pointer<<" "<<
                //    current_sketch.right_child_second_pointer<<endl;
                next_sketch = labels[dir][node_lv][child_first_pointer][child_second_pointer];
                edge_list e( level[node_lv], level[current_end], current_sketch.cost - next_sketch.cost, 
                             current_sketch.weight - next_sketch.weight);
                tmp_path.push_back(e);
                current_end = node_lv;
                current_sketch = next_sketch;
            }
            if(start_v != current_end){
                edge_list e ( level[start_v], level[current_end], current_sketch.cost, current_sketch.weight); 
                tmp_path.push_back(e);
            }
            int l_sz = tmp_path.size();
            for(int k=0; k<tmp_path.size(); k++){
                path.push_back(tmp_path[l_sz -1 -k]);
            }

        }
      

        return;

        int current_end = end_v;
        vector<edge_list>  tmp_path;
        //cout<<"here"<<endl;
        while(current_sketch.pivot_level >=0 ){
            int node_lv = start_v;
            int child_first_pointer, child_second_pointer;
            if(dir==0){
                child_first_pointer = current_sketch.right_child_first_pointer;
                child_second_pointer = current_sketch.right_child_second_pointer;

            }else{
                child_first_pointer = current_sketch.left_child_first_pointer;
                child_second_pointer = current_sketch.left_child_second_pointer;
            }
            //cout<<dir<<endl;
            //cout<<node_lv<<" "<<child_first_pointer<<" "<<child_second_pointer<<" "<<current_sketch.right_child_first_pointer<<" "<<
            //    current_sketch.right_child_second_pointer<<endl;
            next_sketch = labels[dir][node_lv][child_first_pointer][child_second_pointer];
            if(dir == 0){
                edge_list e( level[node_lv], level[current_end], current_sketch.cost - next_sketch.cost, 
                             current_sketch.weight - next_sketch.weight);
                path.push_back(e);
            }else{
                edge_list e(level[current_end], level[node_lv], current_sketch.cost - next_sketch.cost,
                            current_sketch.weight - next_sketch.weight);
                tmp_path.push_back(e);
            }
            triple<int> current_sketch = next_sketch;
            current_end = node_lv;
        }

        if(dir == 1){
            int l_sz = tmp_path.size();
            for(int k=0; k<tmp_path.size(); k++){
                path.push_back(tmp_path[l_sz -1 -k]);
            }
        }
    }


    int c1=0;
    int c2=0;
    
    
    vector<edge_list > constraint_shortest_path_query(int u, int v, int cost){
        vector<edge_list> path;

        if(partition[u] == partition[v]){
            c1++;
            path = dijkstra(u, v, cost);
            return path;
        }else{
            c2++;
            path = cross_partition_query(u,v, cost); 
        }
    }


    void load_order(string order_file_name){
        cout<<"loading order file"<<endl;
        FILE* order_file = fopen(order_file_name.c_str(), "r");
        int vertex;
        
        for(int id=0; fscanf(order_file, "%d", &vertex)==1;id++){
            v2level[vertex]=id;
            level[id] = vertex;
        }
        fclose(order_file);
    }

    void load_index(string index_file_name){
        cout<<"loading index file"<<endl;
        int num_of_labels=0;
        FILE* index_file = fopen(index_file_name.c_str(), "r");
        int id;
        int to_v_id;
        for(int lines =1;fscanf(index_file, "%d: ", &id)==1;lines++){
            int cost, weight, pivot_level, left_child_first_pointer, left_child_second_pointer, right_child_first_pointer, right_child_second_pointer;
            vector<triple<int> > lb_list;
            for(; fscanf(index_file, "%d", &to_v_id)==1; ){
                if(to_v_id == -1){
                    if(lb_list.size()>0){
                        labels[0][id].push_back(lb_list);
                        //label_nodes[0][id].push_back(to_v_id);
                    }
                    break;
                }
                if(fscanf(index_file, " %d %d %d %d %d %d %d", &cost, &weight, &pivot_level,
                          &left_child_first_pointer, &left_child_second_pointer,
                          &right_child_first_pointer, &right_child_second_pointer)!= 7){
                    cout<<"reading index file error in line: "<< lines<<" in reading out-lable"<<endl;
                    exit(0);
                }
                num_of_labels++;
                triple<int> tri;
                tri.v = to_v_id;
                tri.cost = cost;
                tri.weight = weight;
                tri.pivot_level = pivot_level;
                
                tri.left_child_first_pointer = left_child_first_pointer;
                tri.left_child_second_pointer = left_child_second_pointer;
                
                tri.right_child_first_pointer = right_child_first_pointer;
                tri.right_child_second_pointer = right_child_second_pointer;

                
                if(lb_list.size()>0 && lb_list[0].v == to_v_id){
                    lb_list.push_back(tri);
                }else if(lb_list.size()>0){
                    labels[0][id].push_back(lb_list);
                    label_nodes[0][id].push_back(to_v_id);
                    lb_list.clear();
                    lb_list.push_back(tri);
                }else if(lb_list.size() ==0){
                    label_nodes[0][id].push_back(to_v_id);
                    lb_list.push_back(tri);
                }
            }

            lb_list.clear();

            for(; fscanf(index_file, "%d", &to_v_id)==1; ){
                if(to_v_id == -2){
                    if(lb_list.size()>0){
                        labels[1][id].push_back(lb_list);
                    }
                    break;
                }
                if(fscanf(index_file, " %d %d %d %d %d %d %d", &cost, &weight, &pivot_level,
                          &left_child_first_pointer, &left_child_second_pointer,
                          &right_child_first_pointer, &right_child_second_pointer)!= 7){
                    cout<<"reading index file error in line: "<< lines<<" in reading out-lable"<<endl;
                    exit(0);
                }
                num_of_labels++;
                triple<int> tri;
                tri.v = to_v_id;
                tri.cost = cost;
                tri.weight = weight;
                tri.pivot_level = pivot_level;

                tri.left_child_first_pointer = left_child_first_pointer;
                tri.left_child_second_pointer = left_child_second_pointer;
                
                tri.right_child_first_pointer = right_child_first_pointer;
                tri.right_child_second_pointer = right_child_second_pointer;
                if(lb_list.size()>0 && lb_list[0].v == to_v_id){
                    lb_list.push_back(tri);
                }else if(lb_list.size()>0){
                    labels[1][id].push_back(lb_list);
                    label_nodes[1][id].push_back(to_v_id);
                    lb_list.clear();
                    lb_list.push_back(tri);
                }else if(lb_list.size()==0){
                    label_nodes[1][id].push_back(to_v_id);
                    lb_list.push_back(tri);
                }
            }
        }
       
        cout<<"number of labels is: "<<num_of_labels<<endl;
        fclose(index_file);
    }

    int case_1=0;
    int case_1_correct=0;
    int case_2=0;
    int case_2_correct =0;
    pair<int, int> get_left_child_address(triple<int> & father_label, int whose_label, int in_or_out){
        //cout<<"lets' get child address"<<endl;
        int pivot = father_label.pivot_level;
        int leftchild_whose_label, left_child_v_in_label;
        vector<vector<triple<int> > > *checking_label_list_pointer;
        if(in_or_out == 0){
            leftchild_whose_label = pivot <=whose_label? whose_label: pivot;
            left_child_v_in_label = pivot<=whose_label? pivot: whose_label;
            if(pivot<=whose_label){
                checking_label_list_pointer = &labels[0][leftchild_whose_label];
                //case_1++;
            }
            else{
                //case_2++;
                checking_label_list_pointer = &labels[1][leftchild_whose_label];
            }
        }else{
           leftchild_whose_label = pivot;
           left_child_v_in_label = father_label.v;
           checking_label_list_pointer = &labels[1][leftchild_whose_label];
        }
        vector<vector<triple<int> > > & checking_label_list = *checking_label_list_pointer;
        for(int i=0; i<checking_label_list.size(); i++){
            if(checking_label_list[i].size()==0) continue;
            int tmp_v = checking_label_list[i][0].v;
            if(tmp_v == left_child_v_in_label){
                for(int j=0; j<checking_label_list[i].size(); j++){
                    if(checking_label_list[i][j].cost == father_label.left_child_first_pointer &&
                       checking_label_list[i][j].weight == father_label.left_child_second_pointer){
                        //if(pivot<=whose_label) case_1_correct++;
                        //else case_2_correct++;
                        return make_pair(i,j);
                    }
                }
            }
        }
        return make_pair(-1,-1);
        //cout<<"strange"<<endl;
        //assert(0);
    }

    pair<int, int> get_right_child_address(triple<int> & father_label, int whose_label, int in_or_out){
        int pivot = father_label.pivot_level;
        int rightchild_whose_label, right_child_v_in_label;
        vector<vector<triple<int> > > *checking_label_list_pointer;
        if(in_or_out == 0){
            rightchild_whose_label = pivot;
            right_child_v_in_label = father_label.v;
            checking_label_list_pointer = &labels[0][rightchild_whose_label];
        }else{
            rightchild_whose_label = pivot <= whose_label?whose_label: pivot;
            right_child_v_in_label = pivot <= whose_label? pivot : whose_label;
            if(pivot<=whose_label){
                checking_label_list_pointer = &labels[1][rightchild_whose_label];
            }else{
                checking_label_list_pointer = & labels[0][rightchild_whose_label];
            }
        }
        vector<vector<triple<int> > > & checking_label_list = *checking_label_list_pointer;
        for(int i=0; i<checking_label_list.size(); i++){
            if(checking_label_list[i].size()==0) continue;
            int tmp_v = checking_label_list[i][0].v;
            if(tmp_v == right_child_v_in_label){
                for(int j=0; j<checking_label_list[i].size(); j++){
                    if(checking_label_list[i][j].cost == father_label.right_child_first_pointer &&
                       checking_label_list[i][j].weight == father_label.right_child_second_pointer){
                        return make_pair(i,j);
                    }
                }
            }
        }
        return make_pair(-1,-1);
    
    }
    
    
    void children_pointer_linking(){
        int correct=0; int error_num = 0;
        case_1 =0;
        case_1_correct=0;
        case_2 =0;
        case_2_correct=0;
        for(int i=0; i<BG_N; i++){
            for(int j=0; j<labels[0][i].size(); j++){
                for(int k=0; k<labels[0][i][j].size(); k++){
                    int pivot_level = labels[0][i][j][k].pivot_level;
                    if(pivot_level>=0){
                        /*pair<int, int> left_position = get_left_child_address(labels[0][i][j][k], i, 0);
                        if(left_position.first>=0 && left_position.second >=0) case_1_correct++;
                        else{
                            case_1++;
                            cout<<"error when checking "<<i<<"'th label's left children info"<<endl;
                            cout<<"it's vertex id is: "<<level[i]<<endl;
                            cout<<"The label information is as follows"<<endl;
                            triple<int> l_i = labels[0][i][j][k];
                            cout<<"the vertex id in the lable is: "<<level[l_i.v]<<endl;
                            cout<<"the vertex id of pivot is: "<<level[l_i.pivot_level]<<endl;
                            cout<<l_i.v<<" "<<l_i.cost<<" "<<l_i.weight<<" "<<l_i.pivot_level<<" "<<l_i.left_child_first_pointer<<" "<<l_i.left_child_second_pointer<<endl;
                            cout<<l_i.right_child_first_pointer<<" "<<l_i.right_child_second_pointer<<endl;
                            assert(0);
                        }
                        labels[0][i][j][k].left_child_first_pointer = left_position.first;
                        labels[0][i][j][k].left_child_second_pointer = left_position.second;*/
                        pair<int, int> right_position = get_right_child_address(labels[0][i][j][k], i,  0);
                        if(right_position.first>=0 && right_position.second >=0) case_1_correct++;
                        else{
                            case_1++;
                            assert(0);
                        }
                        labels[0][i][j][k].right_child_first_pointer = right_position.first;
                        labels[0][i][j][k].right_child_second_pointer = right_position.second;
                    }
                }
            }
        }

        cout<<"half"<<endl;
        for(int i=0; i<BG_N; i++){
            for(int j=0; j<labels[1][i].size(); j++){
                for(int k=0; k<labels[1][i][j].size(); k++){
                    int pivot_level = labels[1][i][j][k].pivot_level;
                    if(pivot_level>=0){
                        pair<int, int> left_position = get_left_child_address(labels[1][i][j][k], i, 1);
                        if(left_position.first>=0 && left_position.second >=0) case_2_correct++;
                        else{
                            case_2++;
                            assert(0);
                        }
                        labels[1][i][j][k].left_child_first_pointer = left_position.first;
                        labels[1][i][j][k].left_child_second_pointer = left_position.second;
                    }
                }
            }
        }
        cout<<"case 1 error:"<<case_1<<" case_1_correct:"<<case_1_correct<<endl;
        cout<<"case 2 error:"<<case_2<<" case_2_correct:"<<case_2_correct<<endl;
    }

    void output_csp(string filename){
        string full_name = filename+".index"+to_string(APPROX);
        string order_name = filename+".order"+to_string(APPROX);
        FILE* index_file = fopen(full_name.c_str(), "w");
        FILE* order_file = fopen(order_name.c_str(), "w");
        for(int i=0; i<BG_N; i++){
            fprintf(order_file, "%d\n", level[i]);
        }
        cout<<"adding labels"<<endl;
        for(int i=0; i<BG_N; i++){
            if(labels[0][i].size()!=0 || labels[1][i].size()!=0){
                fprintf(index_file, "%d: ", i);
                for(int dir =0; dir<2; dir++){
                    for(unsigned int j=0; j<labels[dir][i].size(); j++){
                        for(unsigned int k=0; k<labels[dir][i][j].size(); k++){
                            fprintf(index_file, "%d %d %d %d %d %d %d %d ", labels[dir][i][j][k].v, labels[dir][i][j][k].cost, labels[dir][i][j][k].weight, 
                                    labels[dir][i][j][k].pivot_level, labels[dir][i][j][k].left_child_first_pointer,labels[dir][i][j][k].left_child_second_pointer, 
                                    labels[dir][i][j][k].right_child_first_pointer, labels[dir][i][j][k].right_child_second_pointer );
                        }
                    }
                    fprintf(index_file, "%d ", (-1-dir) );
                }
                fprintf(index_file, "\n");
            }
        }
        //cout<<"total # of labels: "<<total_labels<<", total # of reduced labels: "<<approximate_reduced_labels<<endl;
        //cout<<"compression ratio: "<<(1.0*(total_labels-approximate_reduced_labels)/(total_labels))<<endl;
        fclose(index_file);
        fclose(order_file);
    }
};

int main(int argc, char** argv){
    string GraphFileName, Graph_Partition_FileName, Graph_Cut_Edges_FileName;
    string Boundary_Graph_FileName;
    string QueryFileName;
  
    string Index_FileName;
    string Order_FileName;

    bool need_boundary_graph = false;
    bool generate_boundary_graph = false;
    bool load_index = false;
    bool do_query_test = false;
    int num_of_queries = 0;
    bool load_order = false; 
    bool output_indices = false;
    cout<<"approximation:"<<APPROX<<endl;
    for ( int i = 1 ; i < argc ; i++ )
	{
		if ( strcmp( argv[i], "-g" ) == 0 ){
			GraphFileName = string(argv[i+1]);
            if(GraphFileName == ""){
                cout<<"Missing GraphFileName"<<endl;
                return 0;
            }
            i++;
		}


		
		if ( strcmp( argv[i], "-gb" ) == 0 ){
            //generate boundary graph
            Graph_Partition_FileName = string(argv[i+1]);
            Graph_Cut_Edges_FileName = string(argv[i+2]);
            need_boundary_graph = true;
            generate_boundary_graph = true;
            i+=2;
		}

        if( strcmp(argv[i], "-lb") == 0 ){
            Boundary_Graph_FileName = string(argv[i+1]);
            Graph_Partition_FileName = string(argv[i+2]);
            Graph_Cut_Edges_FileName = string(argv[i+3]);
            need_boundary_graph = true;
            generate_boundary_graph = false;
            i+=3;
        }

        if(strcmp(argv[i], "-o" ) == 0 ){
            Order_FileName = string(argv[i+1]);
            load_order = true;
            i+=1;
        }

		if ( strcmp( argv[i], "-i" ) == 0 ){
            //Order_FileName = string(argv[i+1]);
            Index_FileName = string(argv[i+1]);
            load_index = true;
            i+=1;
		}
		

		if ( strcmp( argv[i], "-dq" ) == 0 ){
            QueryFileName = string(argv[i+1]);
            num_of_queries = atoi(argv[i+2]);
            do_query_test = 1;
            i+=2;
		}
	}

   CSP csp;

    if(need_boundary_graph == false){
        cout<<"load graph"<<endl;
        csp.load_directed_graph(GraphFileName);
    }else{
        if(generate_boundary_graph == false){
            cout<<"load graph, overlay graph, partition file, and cut edge file"<<endl;
            csp.load_directed_graph(GraphFileName, 
                                    Boundary_Graph_FileName, Graph_Partition_FileName,
                                    Graph_Cut_Edges_FileName);
        }else{
            cout<<"load graph, partition file, and cut edge file"<<endl;
            csp.load_directed_graph(GraphFileName, Graph_Partition_FileName,
                                    Graph_Cut_Edges_FileName);
        }
    }

    if(load_order == true){
        csp.load_order(Order_FileName);
    }else{
        cout<<"ordering vertices"<<endl;
        csp.compute_order_by_coverage_heuristic();
    }

    if(load_index == true){
        //csp.load_order( Order_FileName);
        csp.load_index( Index_FileName); 
    }else{
        cout<<"preprocessing and generating indices"<<endl;
        csp.pruned_Dijkstra(csp.level);
        cout<<"linking children pointer"<<endl;
        csp.children_pointer_linking();
        cout<<"outputing index file"<<endl;
        csp.output_csp(GraphFileName);
    }

   if(do_query_test){
       csp.load_query(QueryFileName);
       cout<<"begin to do query"<<endl;
       clock_t start_q, end_q;
       start_q = clock();
       vector<edge_list> path;
       int reachable_cnt=0;
       for(int i=0; i< num_of_queries ; i++){
           //cout<<"the "<<i<<"'th query: src "<<csp.query_set[i].u<<" dst "<<csp.query_set[i].v<<" cost limit "<<csp.query_set[i].cost_limit <<endl;
           path = csp.constraint_shortest_path_query(csp.query_set[i].u, csp.query_set[i].v, csp.query_set[i].cost_limit);
           if(path.size()>0)
               reachable_cnt++;
       }
       //cout<<csp.c1<<endl;
       //cout<<csp.c2<<endl;
       cout<<"reachable_queries: "<<csp.reachable_queries<<endl;
       cout<<"reachable_cnt: "<<reachable_cnt<<endl;
       cout<<"unfold_shortcut_num: "<<csp.unfold_shortcut_num<<endl;
       cout<<"unfold_edge_num: "<<csp.unfold_edge_num<<endl;
       end_q = clock();
       printf( "Query Time: %0.2lf\n",(double)(end_q-start_q)/CLOCKS_PER_SEC);
   }

}
