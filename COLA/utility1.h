#include <vector>
#include <algorithm>
#include <string>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <queue>
using namespace std;

inline double randrange(double x0, double x1)
{
    return x0 + (x1 - x0) * rand() / ((double) RAND_MAX);
}

struct pWeight
{
	int weight;
	long long p;
	void SetP()
	{
		p = (rand()<<15)+rand();
	}
	bool operator<( const pWeight& w )const
	{
		if ( weight < w.weight ) return true;
		else if ( weight > w.weight ) return false;
		else if ( p == w.p )
		{
			printf("Error in perturbation!!!!!!!\n");
			return true;
		}
		else return ( p < w.p );
	}
	bool operator>( const pWeight& w )const
	{
		if ( weight > w.weight ) return true;
		else if ( weight < w.weight ) return false;
		else if ( p == w.p )
		{
			printf("Error in perturbation!!!!!!!\n");
			return true;
		}
		else return ( p < w.p );
	}
	pWeight operator+( const pWeight& w )
	{
		pWeight tmp;
		tmp.weight = weight + w.weight;
		tmp.p = p + w.p;
		return tmp;
	}
	bool operator==( const pWeight& w )const
	{
		return ( weight==w.weight && p == w.p );
	}
	bool operator==( const int& w )const
	{
		return weight == w;
	}
	bool operator=( const int& w )
	{
		weight = w;
		return true;
	}
	pWeight( int w )
	{
		weight = w;
		p = 0;
	}
	pWeight()
	{
		weight = 0;
		p = 0;
	}
};

struct EMPTYTYPE
{
	bool operator<( const EMPTYTYPE& e )const
	{
		return false;
	}
};

template <typename _T>
struct Triple
{
	int x , y ;
	_T w;
	bool operator<( const Triple& e )const
	{
		if ( x < e.x ) return true;
		else if ( x > e.x ) return false;
		else if ( y < e.y ) return true;
		else if ( y > e.y ) return false;
		else return w < e.w;
	}
	bool operator==( const Triple& e )const
	{
		return ( x == e.x && y == e.y );
	}
	Triple( int xx, int yy, _T ww )
	{
		x=xx;
		y=yy;
		w=ww;
	}
	Triple()
	{
	}
};

template <typename _Key, typename _Value>
struct Key_Value
{
	_Key key;
	_Value value;
	Key_Value( const _Key& k , const _Value& v )
	{
		key = k;
		value = v;
	}
	bool operator==( const Key_Value<_Key,_Value>& p )const
	{
		return ( key == p.key );
	}
	bool operator<( const Key_Value<_Key,_Value>& p )const
	{
		if ( key < p.key ) return true;
		else if ( key > p.key ) return false;
		return ( value < p.value );
	}

	Key_Value( int tmp )
	{
		key = tmp;
		value = -1;
	}

	Key_Value()
	{
	}
};

struct FSCwithCor
{
	int weight, key_point, x, y;
};

struct ZippedFSC
{
	int weight, key_point;
};

struct FinalShortCut
{
	int start, target, key_point, upper, lower;
	int weight;
	bool operator<( const FinalShortCut& e )const
	{
		if ( start < e.start ) return true;
		else if ( start > e.start ) return false;
		else if ( upper > e.upper ) return true;
		else if ( upper < e.upper ) return false;
		else if ( target < e.target ) return true;
		else if ( target > e.target ) return false;
		else if ( weight < e.weight ) return true;
		else if ( weight > e.weight ) return false;
		else return lower < e.lower;
	}
	bool operator==( const FinalShortCut& e )const
	{
		return ( start == e.start && target == e.target );
	}
	Key_Value<int,ZippedFSC> export_JSC()
	{
		ZippedFSC tmp;
		tmp.weight = weight;
		if ( key_point == -1 ) tmp.key_point = (1<<30)+(lower<<25);
		else tmp.key_point = (lower<<25)+key_point;
		return Key_Value<int,ZippedFSC>( (upper<<25)+target , tmp );
	}
	Key_Value<int,FSCwithCor> export_QSC( int x , int y )
	{
		FSCwithCor tmp;
		tmp.weight = weight;
		if ( key_point == -1 ) tmp.key_point = (1<<30);
		else tmp.key_point = key_point;
		tmp.x = x;
		tmp.y = y;
		return Key_Value<int,FSCwithCor>( target , tmp );
	}

	Key_Value<int,ZippedFSC> export_QSC()
	{
		ZippedFSC tmp;
		tmp.weight = weight;
		if ( key_point == -1 ) tmp.key_point = (1<<30);
		else tmp.key_point = key_point;
		return Key_Value<int,ZippedFSC>( target , tmp );
	}

	FinalShortCut()
	{};
	FinalShortCut( int x , int y , int w , int k , int l )
	{
		start = x;
		target = y;
		weight = w;
		key_point = k;
		lower = l;
	}
};

struct WKRR
{
	int weight, keypoint, re[4];
};


struct Level
{
	int* old_level;
	int* new_level;
	int m_NodeNum;
	int a[100][100];

	void free_mem()
	{
		delete[] old_level;
		delete[] new_level;
	}

	Level()
	{
	}
	Level( int n )
	{
		m_NodeNum = n;
		if ( !(old_level = new int[n]) )
		{
			printf("New old level Error!\n");
		}
		if ( !(new_level = new int[n]) )
		{
			printf("New new level Error!\n");
		}
		for ( int i = 0 ; i < n ; ++i )
			old_level[i] = i;
		memset( new_level , 0 , sizeof(int)*n );
	}
	//inline int GetOldLevel( int p )
	//{
	//	return old_level[p];
	//}
	inline int GetNewLevel( int p )
	{
		return new_level[p];
	}
	inline void SetLevel( int p , int l )
	{
		new_level[p] = l;
	}
	//void Refresh()
	//{
	//	memcpy( old_level, new_level, sizeof(int)*m_NodeNum );
	//}

	void OutputDistribution( int r, int range )
	{
		int sum = 0;
		for ( int i = 0 ; i <= range ; ++i )
			for ( int j = 0 ; j <= range ; ++j )
				a[i][j] = 0;
		for ( int i = 0 ; i < m_NodeNum ; ++i )
		{
			a[new_level[i]][old_level[i]]++;
		}
		for ( int i = 0 ; i <= range ; ++i )
		{
			for ( int j = 0 ; j <= range ; ++j )
				if ( a[i][j] > r )
				{
					//printf( "%d %d %d\n", i, j, a[i][j] );
					sum += a[i][j];
				}
		    printf("%d %d\n",sum,m_NodeNum);
			sum = 0;
		}
	}

	void OutputDistribution( int r )
	{
		memset(a , 0 , sizeof(a) );

		for ( int i = 0 ; i < m_NodeNum ; ++i )
			a[new_level[i]][0]++;

		int sum = 0 ;

		for ( int i = 0 ; i <= r ; ++i )
		{
			sum += a[i][0];
			printf("%d %d\n", a[i][0], sum );
		}

	}

	void Output( int L )
	{
		FILE* file = fopen("level.txt","w");
		for ( int i = 0 ; i < m_NodeNum ; ++i )
			if ( new_level[i] < L )
			  fprintf(file,"%d %d %d\n",i,new_level[i],old_level[i]);
		fclose(file);
	}
	void Load( int n, string str )
	{
		new_level = new int[n];
		m_NodeNum = n;
		FILE *file = fopen(str.c_str(),"r");
		int j ,l1,l2;
		for ( int i = 0 ; i < m_NodeNum ; ++i )
		{
			fscanf( file, "%d %d %d", &j , &l1, &l2 );
			if ( i != j ) printf("Not so right!!!\n");
			new_level[j] = l1;
		}
		fclose(file);
	}
	void Load( string str )
	{
		FILE *file = fopen(str.c_str(),"r");
		int j ,l1,l2;
		for ( int i = 0 ; i < m_NodeNum ; ++i )
		{
			fscanf( file, "%d %d %d", &j , &l1, &l2 );
			if ( i != j ) printf("Not so right!!!\n");
			new_level[j] = l1;
			old_level[j] = l2;
		}
		fclose(file);
	}
};

extern const int VectorDefaultSize;

template <typename _T>
class iVector
{
public:
	unsigned int m_size;
	_T* m_data;
	unsigned int m_num;

	void free_mem()
	{
		delete[] m_data;
	}

	iVector()
	{
		//printf("%d\n",VectorDefaultSize);
		m_size = VectorDefaultSize;
		m_data = new _T[VectorDefaultSize];
		m_num = 0;
	}
	iVector( unsigned int n )
	{
		if ( n == 0 )
		{
			n = VectorDefaultSize;
		}
//		printf("iVector allocate: %d\n",n);
		m_size = n;
		m_data = new _T[m_size];
		m_num = 0;
	}

    bool operator==(const iVector& ivec) const{
        if(m_num != ivec.m_num) return false;

        else{
            for(int i=0; i< m_num; i++){
                if(m_data[i] != ivec.m_data[i]) return false;
            }
            return true;
        }
    }

    int size(){
        return m_num;
    }

    void push_back( _T d )
	{
		if ( m_num == m_size )
		{
			re_allocate( m_size*2 );
		}
		m_data[m_num] = d ;
		m_num++;		
	}
	void push_back( const _T* p, unsigned int len )
	{
		while ( m_num + len > m_size )
		{
			re_allocate( m_size*2 );
		}
		memcpy( m_data+m_num, p, sizeof(_T)*len );
		m_num += len;
	}

	void re_allocate( unsigned int size )
	{
		if ( size < m_num )
		{
			return;
		}
		_T* tmp = new _T[size];
		memcpy( tmp, m_data, sizeof(_T)*m_num );
		m_size = size;
		delete[] m_data;
		m_data = tmp;
	}
	void Sort()
	{
		if ( m_num < 20 )
		{
			int k ;
			_T tmp;
			for ( int i = 0 ; i < m_num-1 ; ++i )
			{
				k = i ;
				for ( int j = i+1 ; j < m_num ; ++j )
					if ( m_data[j] < m_data[k] ) k = j ;
				if ( k != i )
				{
					tmp = m_data[i];
					m_data[i] = m_data[k];
					m_data[k] = tmp;
				}
			}
		}
		else sort( m_data, m_data+m_num );
	}
	void unique()
	{
		if ( m_num == 0 ) return;
		Sort();
		unsigned int j = 0;
		for ( unsigned int i = 0 ; i < m_num ; ++i )
			if ( !(m_data[i] == m_data[j]) )
			{
				++j;
				if ( j != i ) m_data[j] = m_data[i];
			}
		m_num = j+1;
	}
	int BinarySearch( _T& data )
	{
		for ( int x = 0 , y = m_num-1 ; x <= y ; )
		{
			int p = (x+y)/2;
			if ( m_data[p] == data ) return p;
			if ( m_data[p] < data ) x = p+1;
			else y = p-1;
		}
		return -1;
	}
	void clean()
	{
		m_num = 0;
	}
	void assign( iVector& t )
	{
		m_num = t.m_num;
		m_size = t.m_size;
		delete[] m_data;
		m_data = t.m_data;
	}

    bool remove_list(int start_idx, int end_idx){
        if(start_idx > end_idx || start_idx<0 || end_idx > m_num){
            printf("wrong parameter, start_idx should be no larger than tend_idx");
            return false;
        }
        int copy_num = m_num - end_idx ;
        memmove(m_data+start_idx, m_data+ end_idx +1, sizeof(_T)*copy_num);
        m_num -= (end_idx - start_idx+1);
        return true;
    }

	bool remove( _T& x )
	{
		for ( int l = 0 , r = m_num ; l < r ; )
		{
			int m = (l+r)/2;

			if ( m_data[m] == x )
			{
                m_num--;
				if ( m_num > m ) memmove( m_data+m, m_data+m+1, sizeof(_T)*(m_num-m) );
				return true;
			}
			else if ( m_data[m] < x ) l = m+1;
			else r = m;
		}
		return false;
	}

	void sorted_insert( _T& x )
	{
		if ( m_num == 0 )
		{
			push_back( x );
			return;
		}

		if ( m_num == m_size ) re_allocate( m_size*2 );

		int l,r;

		for ( l = 0 , r = m_num ; l < r ; )
		{
			int m = (l+r)/2;
			if ( m_data[m] < x ) l = m+1;
			else r = m;
		}

		if ( l < m_num && m_data[l] == x )
		{
            //printf("Insert Duplicate....\n");
            //cout<<x<<endl;
	//		break;
		}
		else
		{
			if ( m_num > l )
			{
				memmove( m_data+l+1, m_data+l, sizeof(_T)*(m_num-l) );
			}
			m_num++;
			m_data[l] = x;
		}
	}

	bool remove_unsorted( _T& x )
	{
		for ( int m = 0 ; m < m_num ; ++m )
		{
			if ( m_data[m] == x )
			{
				m_num--;
				if ( m_num > m ) memcpy( m_data+m, m_data+m+1, sizeof(_T)*(m_num-m) );
				return true;
			}
		}
		return false;
	}

	_T& operator[]( unsigned int i )
	{
		//if ( i < 0 || i >= m_num ) 
		//{
		//	printf("iVector [] out of range!!!\n");
		//}
		return m_data[i];
	}
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//close range check for [] in iVector if release

};

template <typename _T>
struct PairMap
{	
	_T* m_data;
	int m_num;
	int cur;
	iVector<int> occur;
    _T nil;	
    PairMap()
	{
		m_data = NULL;
		m_num = 0;
        //nil = std::make_pair((long)-9,(long)-9);
	    //nil = 1073741834
    }
	void free_mem()
	{
		delete[] m_data;
		occur.free_mem();
	}

	void initialize( int n )
	{
		occur.clean();
		m_num = n;
        nil = std::make_pair((long)-100,(long)-100);
		if ( m_data != NULL )
			delete[] m_data;
		m_data = new _T[m_num];
		for ( int i = 0 ; i < m_num ; ++i )
			m_data[i] = nil;
		cur = 0;
	}
	void clean()
	{
		for ( int i = 0 ; i < occur.m_num ; ++i )
		{
			m_data[occur[i]] = nil;
		}
		occur.clean();
		cur = 0;
	}
	_T get( int p )
	{
		//if ( p < 0 || p >= m_num ) 
		//{
		//	printf("iMap get out of range!!!\n");
		//	return -8;
		//}
		return m_data[p];
	}
	_T& operator[](  int p )
	{
		//if ( i < 0 || i >= m_num ) 
		//{
		//	printf("iVector [] out of range!!!\n");
		//}
		return m_data[p];
	}
	void erase( int p )
	{
		//if ( p < 0 || p >= m_num ) 
		//{
		//	printf("iMap get out of range!!!\n");
		//}
		m_data[p] = nil;
		cur--;
	}
	bool notexist( int p )
	{
		return m_data[p] == nil ;
	}
	bool exist( int p )
	{
        //if(m_data[p] == nil)
        //    printf("fuck\n");
		return !(m_data[p] == nil);
	}
	void insert( int p , _T d )
	{
		//if ( p < 0 || p >= m_num ) 
		//{
		//	printf("iMap insert out of range!!!\n");
		//}
		if ( m_data[p] == nil )
		{
			occur.push_back( p );
			cur++;
		}
        if(d == nil){
            printf("fuck\n");
        }
		m_data[p] = d;
	}
	//close range check when release!!!!!!!!!!!!!!!!!!!!	
};


template <typename _T>
struct iMap
{	
	_T* m_data;
	int m_num;
	int cur;
	iVector<int> occur;
    _T nil;	
    iMap()
	{
		m_data = NULL;
		m_num = 0;
        //nil = std::make_pair((long)-9,(long)-9);
	    //nil = 1073741834;
    }
    iMap(_T tt){
        m_data = NULL;
        m_num = 0;
        nil = tt;
    }
	void free_mem()
	{
		delete[] m_data;
		occur.free_mem();
	}

	void initialize( int n )
	{
		occur.clean();
		m_num = n;
        nil = -9;
		if ( m_data != NULL )
			delete[] m_data;
		m_data = new _T[m_num];
		for ( int i = 0 ; i < m_num ; ++i )
			m_data[i] = nil;
		cur = 0;
	}
	void clean()
	{
		for ( int i = 0 ; i < occur.m_num ; ++i )
		{
			m_data[occur[i]] = nil;
		}
		occur.clean();
		cur = 0;
	}
	_T get( int p )
	{
		//if ( p < 0 || p >= m_num ) 
		//{
		//	printf("iMap get out of range!!!\n");
		//	return -8;
		//}
		return m_data[p];
	}
	_T& operator[](  int p )
	{
		//if ( i < 0 || i >= m_num ) 
		//{
		//	printf("iVector [] out of range!!!\n");
		//}
		return m_data[p];
	}
	void erase( int p )
	{
		//if ( p < 0 || p >= m_num ) 
		//{
		//	printf("iMap get out of range!!!\n");
		//}
		m_data[p] = nil;
		cur--;
	}
	bool notexist( int p )
	{
		return m_data[p] == nil ;
	}
	bool exist( int p )
	{
		return !(m_data[p] == nil);
	}
	void insert( int p , _T d )
	{
		//if ( p < 0 || p >= m_num ) 
		//{
		//	printf("iMap insert out of range!!!\n");
		//}
		if ( m_data[p] == nil )
		{
			occur.push_back( p );
			cur++;
		}
		m_data[p] = d;
	}
	void inc( int p )
	{
		//if ( m_data[p] == nil )
		//{
		//	printf("inc some unexisted point\n");
		//}
		m_data[p]++;
	}
	void inc( int p , int x )
	{
		//if ( m_data[p] == nil )
		//{
		//	printf("inc some unexisted point\n");
		//}
		m_data[p] += x;
	}
	void dec( int p )
	{
		//if ( m_data[p] == nil )
		//{
		//	printf("dec some unexisted point\n" );
		//}
		m_data[p]--;
	}
	//close range check when release!!!!!!!!!!!!!!!!!!!!	
};

struct traversal_pair{
    int v;
    int weight;
    bool operator<(const traversal_pair & tp)  const{
        if(weight > tp.weight) return true;
        else if(weight < tp.weight ) return false;
        else if( v> tp.v) return true;
        else return false;
    }

    bool operator==(const traversal_pair & tp)const{
        return v == tp.v && weight == tp.weight;
    }

    traversal_pair(int vv, int wweight){
        v = vv;
        weight = wweight;
    }
    traversal_pair(){}
};

template <typename __T1>
struct traversal_tuple{
    int v;
    __T1 weight;
    __T1 cost;
    int pivot_level;
    int pivot_pos;
    int to_pivot_weight;
    __T1 to_pivot_cost;
    int hops;
    int direction; //0 denotes path to pivot; 1 denotes path from pivot
    int dominate_paths_minimum_weight;
    bool operator<(const traversal_tuple &tri) const{
        if(cost > tri.cost) return true;
        else if(cost < tri.cost ) return false;
        else if(weight > tri.weight ) return true;
        else if(weight < tri.weight) return false;
        else if ( v > tri.v ) return true;
        else if( v < tri.v ) return false;
        else if( pivot_level > tri.pivot_level) return true;
        else return false;
    }

    bool operator==(const traversal_tuple &tri)const{
        return v == tri.v && weight == tri.weight && 
            cost == tri.cost;
    }

    traversal_tuple(int vv, __T1 wweight, __T1 ccost){
        v = vv;
        weight = wweight;
        cost = ccost;
        pivot_level = -2;
        direction = -1;
        hops = 0;
        dominate_paths_minimum_weight = -9;
    }
    traversal_tuple(){}
};


template<typename _T1, typename _T2>
struct pareto_pair{
    _T1 first;
    _T2 second;
    int pivot_level;
    int pivot_pos;
    _T1 pivot_first;
    _T2 pivot_second;
    int hops;
    int dominate_paths_minimum_weight;
    int expanded;
    int direction;//0 denotes paths to pivot; 1 denotes path from pivot
    pareto_pair(_T1 cost, _T2 weight){
        first = cost;
        second = weight;
        pivot_level = -2;
        //is_domiated = false;
        pivot_first = -1;
        pivot_second = -1;
        direction = -1;
        hops = 0;
        expanded = 0;
        dominate_paths_minimum_weight = -9;
    }

    bool operator<(const pareto_pair pp) const{
        if(first <pp.first) return true;
        else if(first > pp.first) return false;
        else if(second < pp.second) return true;
        else if(second > pp.second) return false;
        else if( pivot_level < pp.pivot_level) return true;
        else return false;
    }

    bool operator==(const pareto_pair pp)const{
        return first == pp.first && second == pp.second;
    }
};

template <typename _T>
struct MMap
{	
	_T* m_data;
	int m_num;
	int cur;
	iVector<int> occur;
    _T nil;	
    MMap()
	{
		m_data = NULL;
		m_num = 0;
    }
    MMap(_T tt){
        m_data = NULL;
        m_num = 0;
        nil = tt;
    }
	void free_mem()
	{
		delete[] m_data;
		occur.free_mem();
	}

    void set_nil(_T& nill){
        nil = nill;
    }

	void initialize( int n )
	{
		occur.clean();
		m_num = n;
        //vector<pair<int, int> > tmp;
        //tmp.push_back(make_pair(-1,-1)); 
        //iVector<pair<int, int> > tmp;
        //tmp.push_back(make_pair(-1,-1) );
        //nil = tmp;
		if ( m_data != NULL )
			delete[] m_data;
		m_data = new _T[m_num];
		for ( int i = 0 ; i < m_num ; ++i )
			m_data[i] = nil;
		cur = 0;
	}
	void clean()
	{
		for ( int i = 0 ; i < occur.m_num ; ++i )
		{
			m_data[occur[i]] = nil;
		}
		occur.clean();
		cur = 0;
	}
	_T get( int p )
	{
		return m_data[p];
	}
	_T& operator[](  int p )
	{
		return m_data[p];
	}
	void erase( int p )
	{
		m_data[p] = nil;
		cur--;
	}
	bool notexist( int p )
	{
		return m_data[p] == nil ;
	}
	bool exist( int p )
	{
		return !(m_data[p] == nil);
	}
	void insert( int p , _T d )
	{
		if ( m_data[p] == nil )
		{
			occur.push_back( p );
			cur++;
		}
		m_data[p] = d;
	}
	void inc( int p )
	{
		m_data[p]++;
	}
	void inc( int p , int x )
	{
		m_data[p] += x;
	}
	void dec( int p )
	{
		m_data[p]--;
	}
	//close range check when release!!!!!!!!!!!!!!!!!!!!	
};

template <typename _T>
struct Comparable_Map
{	
	_T* m_data;
	int m_num;
	int cur;
	iVector<int> occur;
    _T nil;	
    Comparable_Map()
	{
		m_data = NULL;
		m_num = 0;
    }
    Comparable_Map(_T tt){
        m_data = NULL;
        m_num = 0;
        nil = tt;
    }
	void free_mem()
	{
		delete[] m_data;
		occur.free_mem();
	}

    void set_nil(_T& nill){
        nil = nill;
    }

	void initialize( int n )
	{
		occur.clean();
		m_num = n;
        //vector<pair<int, int> > tmp;
        //tmp.push_back(make_pair(-1,-1)); 
        //iVector<pair<int, int> > tmp;
        //tmp.push_back(make_pair(-1,-1) );
        //nil = tmp;
		if ( m_data != NULL )
			delete[] m_data;
		m_data = new _T[m_num];
		for ( int i = 0 ; i < m_num ; ++i )
			m_data[i] = nil;
		cur = 0;
	}
	void clean()
	{
		for ( int i = 0 ; i < occur.m_num ; ++i )
		{
			m_data[occur[i]] = nil;
		}
		occur.clean();
		cur = 0;
	}
	_T get( int p )
	{
		return m_data[p];
	}
	_T& operator[](  int p )
	{
		return m_data[p];
	}
	void erase( int p )
	{
		m_data[p] = nil;
		cur--;
	}
	bool notexist( int p )
	{
		return m_data[p] == nil ;
	}
	bool exist( int p )
	{
		return !(m_data[p]==nil);
	}
	void insert( int p , _T d )
	{
		if ( m_data[p] == nil )
		{
			occur.push_back( p );
			cur++;
		}
		m_data[p] = d;
	}
	void inc( int p )
	{
		m_data[p]++;
	}
	void inc( int p , int x )
	{
		m_data[p] += x;
	}
	void dec( int p )
	{
		m_data[p]--;
	}
	//close range check when release!!!!!!!!!!!!!!!!!!!!	
};
template <typename _T>
struct VecMap
{	
	_T* m_data;
	int m_num;
	int cur;
	iVector<int> occur;
    _T nil;	
    VecMap()
	{
		m_data = NULL;
		m_num = 0;
        //nil = std::make_pair((long)-9,(long)-9);
	    //nil = 1073741834;
    }
    VecMap(_T tt){
        m_data = NULL;
        m_num = 0;
        nil = tt;
    }
	void free_mem()
	{
		delete[] m_data;
		occur.free_mem();
	}

	void initialize( int n )
	{
		occur.clean();
		m_num = n;
        vector<pareto_pair<int,int> > tmp;
        pareto_pair<int,int> pp(-1,-1);
        tmp.push_back(pp);
        //vector<pair<int, int> > tmp;
        //tmp.push_back(make_pair(-1,-1)); 
        //iVector<pair<int, int> > tmp;
        //tmp.push_back(make_pair(-1,-1) );
        nil = tmp;
		if ( m_data != NULL )
			delete[] m_data;
		m_data = new _T[m_num];
		for ( int i = 0 ; i < m_num ; ++i )
			m_data[i] = nil;
		cur = 0;
	}
	void clean()
	{
		for ( int i = 0 ; i < occur.m_num ; ++i )
		{
			m_data[occur[i]] = nil;
		}
		occur.clean();
		cur = 0;
	}
	_T get( int p )
	{
		//if ( p < 0 || p >= m_num ) 
		//{
		//	printf("iMap get out of range!!!\n");
		//	return -8;
		//}
		return m_data[p];
	}
	_T& operator[](  int p )
	{
		//if ( i < 0 || i >= m_num ) 
		//{
		//	printf("iVector [] out of range!!!\n");
		//}
		return m_data[p];
	}
	void erase( int p )
	{
		//if ( p < 0 || p >= m_num ) 
		//{
		//	printf("iMap get out of range!!!\n");
		//}
		m_data[p] = nil;
		cur--;
	}
	bool notexist( int p )
	{
		return m_data[p] == nil ;
	}
	bool exist( int p )
	{
		return !(m_data[p] == nil);
	}
	void insert( int p , _T d )
	{
		//if ( p < 0 || p >= m_num ) 
		//{
		//	printf("iMap insert out of range!!!\n");
		//}
		if ( m_data[p] == nil )
		{
			occur.push_back( p );
			cur++;
		}
		m_data[p] = d;
	}
	void inc( int p )
	{
		//if ( m_data[p] == nil )
		//{
		//	printf("inc some unexisted point\n");
		//}
		m_data[p]++;
	}
	void inc( int p , int x )
	{
		//if ( m_data[p] == nil )
		//{
		//	printf("inc some unexisted point\n");
		//}
		m_data[p] += x;
	}
	void dec( int p )
	{
		//if ( m_data[p] == nil )
		//{
		//	printf("dec some unexisted point\n" );
		//}
		m_data[p]--;
	}
	//close range check when release!!!!!!!!!!!!!!!!!!!!	
};
struct dMap
{	
	double* m_data;
	int m_num;
	int cur;
	iVector<int> occur;
	void initialize( int n )
	{
		occur.re_allocate(1000);
		m_num = n;
		m_data = new double[m_num];
		for ( int i = 0 ; i < m_num ; ++i )
			m_data[i] = 1e50;
		cur = 0;
	}
	void clean()
	{
		for ( unsigned int i = 0 ; i < occur.m_num ; ++i )
		{
			m_data[occur[i]] = 1e50;
		}
		occur.clean();
		cur = 0;
	}
	double get( int p )
	{
		//if ( p < 0 || p >= m_num ) 
		//{
		//	printf("iMap get out of range!!!\n");
		//	return -8;
		//}
		return m_data[p];
	}
	void erase( int p )
	{
		//if ( p < 0 || p >= m_num ) 
		//{
		//	printf("iMap get out of range!!!\n");
		//}
		m_data[p] = 1e50;
		cur--;
	}
	bool notexist( int p )
	{
		return m_data[p] == 1e50 ;
	}
	bool exist( int p )
	{
		return m_data[p] != 1e50;
	}
	void insert( int p , int d )
	{
		//if ( p < 0 || p >= m_num ) 
		//{
		//	printf("iMap insert out of range!!!\n");
		//}
		if ( m_data[p] == 1e50 )
		{
			occur.push_back( p );
			cur++;
		}
		m_data[p] = d;
	}
	//close range check when release!!!!!!!!!!!!!!!!!!!!	
};

template<typename _Value>
struct Pareto_Queue
{
    std::priority_queue<_Value> pq;
    Pareto_Queue(){}
    Pareto_Queue(priority_queue<_Value> pqq){
        pq = pqq;
    }
    bool operator<(const priority_queue<_Value>& pq2){
        if(pq.size()==0) return false;
        else if(pq.size()==0) return true;
        else if(pq.top() <pq2.top()) return true;
        else return false;
    }
};

template <typename _Value>
struct Pareto_Heap
{
	iMap<int> pos;
	iVector<Key_Value<int,priority_queue<_Value> > > m_data;

	void initialize( int n )
	{
		pos.initialize(n);		
		m_data.re_allocate(1000);
	}

	int head()
	{
		return m_data[0].key.top();
	}

	void clean()
	{
		pos.occur.clean();
	}
    
    void free_mem()
    {
        pos.free_mem();
        m_data.free_mem();
    }
	void DeepClean()
	{
		pos.clean();
		m_data.clean();
	}

	void insert( int x, _Value y )
	{
		if ( pos.notexist(x) )
		{
            priority_queue<_Value> pq;
            pq.push(y);
			m_data.push_back( Key_Value<int,_Value>(x,pq) );
            pos.insert( x , m_data.m_num-1 );
			up( m_data.m_num-1 );
		}
		else
		{
           m_data[pos.get(x)].value.push_back(y);
           if(m_data[pos.get(x)].value.top() == y)
               up(pos.get(x));

			/*if ( y < m_data[pos.get(x)].value )
			{
				m_data[ pos.get(x) ].value = y;
				up( pos.get(x) );
			}
			else
			{
				m_data[ pos.get(x) ].value = y;
				down( pos.get(x) );
			}*/
		}
	}
	_Value pop()
	{
        priority_queue<_Value>& pq = m_data[0].value;
        _Value ret_val = pq.top();
        pq.pop();
        if(!pq.empty()){
            _Value remove_p = pq.top();
            while(!pq.empty() && remove_p.cost >=ret_val.cost && remove_p.weight >= ret_val.weight){
                pq.pop();
                remove_p = pq.top();
            }
        }

        if(pq.size()==0){
            int tmp = m_data[0].key;
            pos.erase(tmp);
            m_data.m_num--;
            return ret_val;
        }
		down(0);
		return ret_val;
	}
	bool empty()
	{
		return ( m_data.m_num == 0 );
	}
	void up( int p )
	{
		Key_Value<int,_Value> x = m_data[p];
		for ( ; p > 0 && x.value.top() < m_data[(p-1)/2].value.top() ; p = (p-1)/2 )
		{
			m_data[p] = m_data[(p-1)/2];
			pos.insert(m_data[p].key, p);
		}
		m_data[p] = x;
		pos.insert(x.key,p);
	}
	void down( int p )
	{
		//Key_Value<int,int> x = m_data[m_data.m_num-1];
		//m_data.m_num--;

		Key_Value<int,_Value> tmp;

		for ( int i ; p < m_data.m_num ; p = i )
		{
			//m_data[p] = x;
			if ( p*2+1 < m_data.m_num && m_data[p*2+1].value.top() < m_data[p].value.top() )
				i = p*2+1;
			else i = p;
			if ( p*2+2 < m_data.m_num && m_data[p*2+2].value.top() < m_data[i].value.top() )
				i = p*2+2;
			if ( i == p ) break;

			tmp = m_data[p];
			m_data[p] = m_data[i];
			pos.insert( m_data[p].key, p );
			m_data[i] = tmp;
		}

		pos.insert(m_data[p].key,p);
	}

};
template <typename _Value>
struct iHeap
{
	iMap<int> pos;
	iVector<Key_Value<int,_Value> > m_data;

	void initialize( int n )
	{
		pos.initialize(n);		
		m_data.re_allocate(1000);
	}

	int head()
	{
		return m_data[0].key;
	}

	void clean()
	{
		pos.occur.clean();
	}
    
    void free_mem()
    {
        pos.free_mem();
        m_data.free_mem();
    }
	void DeepClean()
	{
		pos.clean();
		m_data.clean();
	}

	void insert( int x, _Value y )
	{
		if ( pos.notexist(x) )
		{
			m_data.push_back( Key_Value<int,_Value>(x,y) );
			pos.insert( x , m_data.m_num-1 );
			up( m_data.m_num-1 );
		}
		else
		{
			if ( y < m_data[pos.get(x)].value )
			{
				m_data[ pos.get(x) ].value = y;
				up( pos.get(x) );
			}
			else
			{
				m_data[ pos.get(x) ].value = y;
				down( pos.get(x) );
			}
		}
	}
	int pop()
	{
		int tmp = m_data[0].key;
		pos.erase(tmp);
		if ( m_data.m_num == 1 ) 
		{
			m_data.m_num--;
			return tmp;
		}
		m_data[0] = m_data[m_data.m_num-1];
		m_data.m_num--;
		down(0);
		return tmp;
	}
	bool empty()
	{
		return ( m_data.m_num == 0 );
	}
	void up( int p )
	{
		Key_Value<int,_Value> x = m_data[p];
		for ( ; p > 0 && x.value < m_data[(p-1)/2].value ; p = (p-1)/2 )
		{
			m_data[p] = m_data[(p-1)/2];
			pos.insert(m_data[p].key, p);
		}
		m_data[p] = x;
		pos.insert(x.key,p);
	}
	void down( int p )
	{
		//Key_Value<int,int> x = m_data[m_data.m_num-1];
		//m_data.m_num--;

		Key_Value<int,_Value> tmp;

		for ( int i ; p < m_data.m_num ; p = i )
		{
			//m_data[p] = x;
			if ( p*2+1 < m_data.m_num && m_data[p*2+1].value < m_data[p].value )
				i = p*2+1;
			else i = p;
			if ( p*2+2 < m_data.m_num && m_data[p*2+2].value < m_data[i].value )
				i = p*2+2;
			if ( i == p ) break;

			tmp = m_data[p];
			m_data[p] = m_data[i];
			pos.insert( m_data[p].key, p );
			m_data[i] = tmp;
		}

		pos.insert(m_data[p].key,p);
	}

};

struct Traversal_Heap
{
	iMap<int> pos;
	iVector<Key_Value<int,traversal_tuple<int>> > m_data;

	void initialize( int n )
	{
		pos.initialize(n);		
		m_data.re_allocate(1000);
	}

	int head()
	{
		return m_data[0].key;
	}

    traversal_tuple<int> top()
    {
        return m_data[0].value;
    }

	void clean()
	{
		pos.occur.clean();
	}
    
    void free_mem()
    {
        pos.free_mem();
        m_data.free_mem();
    }
	void DeepClean()
	{
		pos.clean();
		m_data.clean();
	}

    bool smaller_than(const traversal_tuple<int>& t1, const traversal_tuple<int>& t2){
        if(t1.cost < t2.cost) return true;
        else if(t1.cost > t2.cost) return false;
        else if(t1.weight < t2.weight) return true;
        else if(t1.weight > t2.weight) return false;
        else if(t1.v < t2.v ) return true;
        else return false;
    }

	void insert( int x, traversal_tuple<int> y )
	{
		if ( pos.notexist(x) )
		{
			m_data.push_back( Key_Value<int,traversal_tuple<int> >(x,y) );
			pos.insert( x , m_data.m_num-1 );
			up( m_data.m_num-1 );
		}
		else
		{
			if ( smaller_than( y,  m_data[pos.get(x)].value)  )
			{
				m_data[ pos.get(x) ].value = y;
				up( pos.get(x) );
			}
			else
			{
				m_data[ pos.get(x) ].value = y;
				down( pos.get(x) );
			}
		}
	}
	int pop()
	{
		int tmp = m_data[0].key;
		pos.erase(tmp);
		if ( m_data.m_num == 1 ) 
		{
			m_data.m_num--;
			return tmp;
		}
		m_data[0] = m_data[m_data.m_num-1];
		m_data.m_num--;
		down(0);
		return tmp;
	}
	bool empty()
	{
		return ( m_data.m_num == 0 );
	}
	void up( int p )
	{
		Key_Value<int,traversal_tuple<int> > x = m_data[p];
		for ( ; p > 0 && smaller_than(x.value , m_data[(p-1)/2].value) ; p = (p-1)/2 )
		{
			m_data[p] = m_data[(p-1)/2];
			pos.insert(m_data[p].key, p);
		}
		m_data[p] = x;
		pos.insert(x.key,p);
	}
	void down( int p )
	{
		//Key_Value<int,int> x = m_data[m_data.m_num-1];
		//m_data.m_num--;

		Key_Value<int, traversal_tuple<int> > tmp;

		for ( int i ; p < m_data.m_num ; p = i )
		{
			//m_data[p] = x;
			if ( p*2+1 < m_data.m_num && smaller_than( m_data[p*2+1].value , m_data[p].value) )
				i = p*2+1;
			else i = p;
			if ( p*2+2 < m_data.m_num && smaller_than ( m_data[p*2+2].value , m_data[i].value) )
				i = p*2+2;
			if ( i == p ) break;

			tmp = m_data[p];
			m_data[p] = m_data[i];
			pos.insert( m_data[p].key, p );
			m_data[i] = tmp;
		}

		pos.insert(m_data[p].key,p);
	}

};
/*
template<typename T1, typename T2>
struct SkyLineHeap
{
	iMap<int> pos;
	iVector<Key_Value<int, vector<pair<T1, T2> > > > m_data;

	void initialize( int n )
	{
		pos.initialize(n);		
		m_data.re_allocate(1000);
	}

	int head()
	{
		return m_data[0].key;
	}

	void clean()
	{
		pos.occur.clean();
	}
    
    void free_mem()
    {
        pos.free_mem();
        m_data.free_mem();
    }
	void DeepClean()
	{
		pos.clean();
		m_data.clean();
	}

	void insert( int x, pair<T1, T2> y )
	{
		if ( pos.notexist(x) )
		{
            vector<pair< T1, T2> > list;
            list.push_back(y);
			m_data.push_back( Key_Value<int,vector<pair<T1,T2> > >(x,list) );
			pos.insert( x , m_data.m_num-1 );
			up( m_data.m_num-1 );
		}
		else
		{
            vector<pair<T1, T2>> &x_list = m_data[pos.get(x)].value;
            int l, r;
		    int start_idx = -1;
            int end_idx = x_list.m_num;
            int dominated = false;
            int need_update = false;
            int insert_pos=-1;
            for(l=0,r=x_list.size(); l<r;){
                int m = l+ (r-l)/2;
                if(x_list[m] < p) l=m+1;
                else r = m;
            }
            insert_pos = l;
            if(insert_pos>0){
                if(x_list[insert_pos-1].first <= y.first 
                   && x_list[insert_pos-1].second <=y.second)
                    dominated =true;
            }
            if(dominated) return;

            for(int w=l; w<x_list.m_num; w++){
                if(x_list[w].first>=y.first && x_list[w].second >= y.second){
                    if(need_update == false){
                        need_update = true;
                        start_idx = w;
                    }
                    end_idx = w;
                }else if(need_update == true) break;
            }
            
            if(!need_update){
                typename vector<pair<T1,T2> >::iterator it = x_list.begin()+l;
                x_list[w].insert(it, y);
            }else{
                assert(l == start_idx);
                if(start_idx == end_idx)
                    x_list[l] = y;
                else{
                    assert(start_idx < end_idx);
                    x_list[l] = y;
                    vector<pair<T1, T2> >::iterator it = x_list.begin() + l;
                    x_list.erase( it+1, it+end_idx);
                }
            }

            if(insert_pos ==0){
               up(pos.get(x)); 
            }
		}
	}
	int pop()
	{
		int tmp = m_data[0].key;
		pos.erase(tmp);
		if ( m_data.m_num == 1 ) 
		{
			m_data.m_num--;
			return tmp;
		}
		m_data[0] = m_data[m_data.m_num-1];
		m_data.m_num--;
		down(0);
		return tmp;
	}
	bool empty()
	{
		return ( m_data.m_num == 0 );
	}
	void up( int p )
	{
		Key_Value<int,_Value> x = m_data[p];
		for ( ; p > 0 && x.value < m_data[(p-1)/2].value ; p = (p-1)/2 )
		{
			m_data[p] = m_data[(p-1)/2];
			pos.insert(m_data[p].key, p);
		}
		m_data[p] = x;
		pos.insert(x.key,p);
	}
	void down( int p )
	{
		//Key_Value<int,int> x = m_data[m_data.m_num-1];
		//m_data.m_num--;

		Key_Value<int,_Value> tmp;

		for ( int i ; p < m_data.m_num ; p = i )
		{
			//m_data[p] = x;
			if ( p*2+1 < m_data.m_num && m_data[p*2+1].value < m_data[p].value )
				i = p*2+1;
			else i = p;
			if ( p*2+2 < m_data.m_num && m_data[p*2+2].value < m_data[i].value )
				i = p*2+2;
			if ( i == p ) break;

			tmp = m_data[p];
			m_data[p] = m_data[i];
			pos.insert( m_data[p].key, p );
			m_data[i] = tmp;
		}

		pos.insert(m_data[p].key,p);				
	}

};
*/

struct kvHash
{
	int m_size;
	Key_Value<int,int>* m_data;

	kvHash( int n )
	{
		m_size = n;
		m_data = new Key_Value<int,int>[m_size];
	}
	
	kvHash()
	{
	}

	void free_mem()
	{
		delete[] m_data;
	}

	void initialize( int n )
	{
		m_size = n;
		m_data = new Key_Value<int,int>[m_size];
		clean();
	}

	int find( int x )
	{
		int p = x % m_size;

		while ( m_data[p].key > -1 )
		{
			if ( m_data[p].key == x ) return p;
			p++;
			if ( p == m_size ) p = 0;
		}

		return p;
	}

	int get( int x )
	{
		int p = find( x );
		if ( m_data[p].key != x )
			printf("Error in getting non-exist elements!!!!\n");
		return m_data[p].value;
	}

	bool notexist( int x )
	{
		return !exist(x);
	}

	bool exist( int x )
	{
		return (m_data[find(x)].key == x);
	}

	void insert( int x , int y )
	{
		int p = find(x);
		m_data[p].key = x;
		m_data[p].value = y;
	}

	void erase( int x )
	{
		m_data[find(x)].key = -1;
	}

	void clean()
	{
		for ( int i = 0 ; i < m_size ; ++i )
			m_data[i].key = -1;
	}
};


template<typename _Value>
struct hHeap
{
	kvHash pos;
	iVector<Key_Value<int,_Value> > m_data;

	void initialize( int n )
	{
		pos.initialize(n);		
	}

	void free_mem()
	{
		pos.free_mem();
		m_data.free_mem();
	}

	int head()
	{
		return m_data[0].key;
	}

	void DeepClean()
	{
		pos.clean();
		m_data.clean();
	}

	void insert( int x, _Value y )
	{
		if ( pos.notexist(x) )
		{
			m_data.push_back( Key_Value<int,_Value>(x,y) );
			pos.insert( x , m_data.m_num-1 );
			up( m_data.m_num-1 );
		}
		else
		{
			if ( y < m_data[pos.get(x)].value )
			{
				m_data[ pos.get(x) ].value = y;
				up( pos.get(x) );
			}
			else
			{
				m_data[ pos.get(x) ].value = y;
				down( pos.get(x) );
			}
		}
	}
	int pop()
	{
		int tmp = m_data[0].key;
		pos.erase(tmp);
		if ( m_data.m_num == 1 ) 
		{
			m_data.m_num--;
			return tmp;
		}
		m_data[0] = m_data[m_data.m_num-1];
		m_data.m_num--;
		down(0);
		return tmp;
	}
	bool empty()
	{
		return ( m_data.m_num == 0 );
	}
	void up( int p )
	{
		Key_Value<int,_Value> x = m_data[p];
		for ( ; p > 0 && x.value < m_data[(p-1)/2].value ; p = (p-1)/2 )
		{
			m_data[p] = m_data[(p-1)/2];
			pos.insert(m_data[p].key, p);
		}
		m_data[p] = x;
		pos.insert(x.key,p);
	}
	void down( int p )
	{
		//Key_Value<int,int> x = m_data[m_data.m_num-1];
		//m_data.m_num--;

		Key_Value<int,_Value> tmp;

		for ( int i ; p < m_data.m_num ; p = i )
		{
			//m_data[p] = x;
			if ( p*2+1 < m_data.m_num && m_data[p*2+1].value < m_data[p].value )
				i = p*2+1;
			else i = p;
			if ( p*2+2 < m_data.m_num && m_data[p*2+2].value < m_data[i].value )
				i = p*2+2;
			if ( i == p ) break;

			tmp = m_data[p];
			m_data[p] = m_data[i];
			pos.insert( m_data[p].key, p );
			m_data[i] = tmp;
		}

		pos.insert(m_data[p].key,p);				
	}

};

template<typename _Value>
struct bHeap
{
	int* pos;
	iVector<Key_Value<int,_Value> > m_data;
	int base;

	void initialize( int n , int b )
	{
		base = b;
		pos = new int[n];
		for ( int i = 0 ; i < n ; ++i )
			pos[i] = -1;
	}

	void free_mem()
	{
		delete[] pos;
		m_data.free_mem();
	}

	int head()
	{
		return m_data[0].key;
	}

	void DeepClean()
	{
		for ( int i = 0 ; i < m_data.m_num ; ++i )
			pos[m_data[i].key-base] = -1;
		m_data.clean();
	}

	void insert( int x, _Value y )
	{
		if ( pos[x-base] == -1 )
		{
			m_data.push_back( Key_Value<int,_Value>(x,y) );
			pos[x-base] = m_data.m_num-1;
			up( m_data.m_num-1 );
		}
		else
		{
			if ( y < m_data[pos[x-base]].value )
			{
				m_data[ pos[x-base] ].value = y;
				up( pos[x-base] );
			}
			else
			{
				m_data[ pos[x-base] ].value = y;
				down( pos[x-base] );
			}
		}
	}
	int pop()
	{
		int tmp = m_data[0].key;
		pos[tmp-base] = -1;
		if ( m_data.m_num == 1 ) 
		{
			m_data.m_num--;
			return tmp;
		}
		m_data[0] = m_data[m_data.m_num-1];
		m_data.m_num--;
		down(0);
		return tmp;
	}
	bool empty()
	{
		return ( m_data.m_num == 0 );
	}
	void up( int p )
	{
		Key_Value<int,_Value> x = m_data[p];
		for ( ; p > 0 && x.value < m_data[(p-1)/2].value ; p = (p-1)/2 )
		{
			m_data[p] = m_data[(p-1)/2];
			pos[m_data[p].key-base] = p;
		}
		m_data[p] = x;
		pos[x.key-base] = p;
	}
	void down( int p )
	{
		//Key_Value<int,int> x = m_data[m_data.m_num-1];
		//m_data.m_num--;

		Key_Value<int,_Value> tmp;

		for ( int i ; p < m_data.m_num ; p = i )
		{
			//m_data[p] = x;
			if ( p*2+1 < m_data.m_num && m_data[p*2+1].value < m_data[p].value )
				i = p*2+1;
			else i = p;
			if ( p*2+2 < m_data.m_num && m_data[p*2+2].value < m_data[i].value )
				i = p*2+2;
			if ( i == p ) break;

			tmp = m_data[p];
			m_data[p] = m_data[i];
			pos[m_data[p].key-base] = p;
			m_data[i] = tmp;
		}

		pos[m_data[p].key-base] = p;				
	}

};



struct rHeap
{
	iVector<Key_Value<int,int> > m_data;
	void insert( int x, int y )
	{
		m_data.push_back( Key_Value<int,int>(x,y) );
		up( m_data.m_num-1 );
	}
	int pop()
	{
		if ( m_data.m_num == 0 ) 
		{
			printf( "Heap empty when pop!!!!\n" );
			return 1073741834;
		}
		int tmp = m_data[0].key;
		down();
		return tmp;
	}
	bool empty()
	{
		return ( m_data.m_num == 0 );
	}
	void up( int p )
	{
		Key_Value<int,int> x = m_data[p];
		for ( ; p > 0 && x.value < m_data[(p-1)/2].value ; p = (p-1)/2 )
		{
			m_data[p] = m_data[(p-1)/2];
		}
		m_data[p] = x;
	}
	void down()
	{
		if ( m_data.m_num == 1 ) 
		{
			m_data.m_num--;
			return;
		}
		Key_Value<int,int> x = m_data[m_data.m_num-1];
		m_data.m_num--;
		int p = 0, i;

		for ( ; p < m_data.m_num ; p = i )
		{
			m_data[p] = x;
			if ( p*2+1 < m_data.m_num && m_data[p*2+1].value < m_data[p].value )
				i = p*2+1;
			else i = p;
			if ( p*2+2 < m_data.m_num && m_data[p*2+2].value < m_data[i].value )
				i = p*2+2;
			if ( i == p ) break;
			m_data[p] = m_data[i];
		}

	}
};


class DoubleHash
{
	int bset;
	int R, U;
	Triple<int>* m_data;
public:	
	void clean()
	{
		for ( int i = 0 ; i < R ; ++i )
			m_data[i].w = -1;
		U = 0;
	}

	DoubleHash( int range )
	{
		int size = 1 ;
		for ( bset = 0 ; size < range ; size *= 2, bset++ );

		bset = 31-bset;

		R = 128*size-1;
		m_data = new Triple<int>[R];

		clean();
	}

	int f( int x , int y )
	{
		return ((x<<bset)|y)%R;
	}

	int find( int x , int y )
	{
		int i;
		for ( i = f(x,y); m_data[i].w > -1 ; ++i )
		{			
			if ( x == m_data[i].x && y == m_data[i].y ) return m_data[i].w;
			if ( i == R-1 ) i = -1;
		}
		return -1;
	}

	void insert( int x , int y , int w )
	{
		U++;
		int i;
		for ( i = f(x,y); m_data[i].w > -1 ; ++i )
		{
			if ( i == R-1 ) i = -1;
		}
		m_data[i].w = w ;
		m_data[i].x = x ;
		m_data[i].y = y ;
	}

	void usage()
	{
		printf("Hashing usage: %0.2lf\n", (double)U/(double)R);
	}
};

struct PendingQueue
{
	iVector<int> queue;
	iMap<int> pos;
	int point;

	PendingQueue( int n )
	{
		point = 0;
		pos.initialize( n );
	}

	void clean()
	{
		queue.clean();
		point = 0;
		pos.clean();
	}

	void push_back( int x )
	{
		pos.insert( x , queue.m_num );
		queue.push_back( x );
	}

	int pop()
	{
		for ( ; point < queue.m_num ; point++ )
		{
			if ( pos.get(queue[point]) == point )
			{
				point++;
				return queue[point-1];
			}
		}
		return -1;
	}
};

extern int CodeEdgeRange;

struct CodeEdge
{
	int x , y , jump;
	pWeight weight;
	unsigned short code;
	CodeEdge()
	{
		x = CodeEdgeRange;
		y = 0;		
		code = 0;
	}
	bool CodeCompare( const unsigned short c1 , const unsigned short c2 )const
	{
		int r = 0;
		unsigned short p = 1;
		for ( int i = 0 ; i < 16 ; ++i )
		{
			if ( c1&p ) r++;
			if ( c2&p ) r--;
			if ( i < 15 ) p *= 2;
		}
		return r >= 0;
	}
	bool operator<( const CodeEdge& e )const
	{
		if ( x < e.x ) return true;
		else if ( x > e.x ) return false;
		else if ( y < e.y ) return true;
		else if ( y > e.y ) return false;
		else return code < e.code;
		//else return CodeCompare( code , e.code );
	}
	bool operator==( const CodeEdge& e )const
	{
		return ( x == e.x && y == e.y );
	}
};

extern const int HashDefaultRelativeSize;

struct iEdgeHash
{
	int *m_start;
	iVector<CodeEdge> hash_table;
	iVector<int> hash_map;
	int m_num;
	int m_EdgeNum;
	int occupied;
	int AverageSpace;

	int find( int x , int y )
	{
		for ( int i = GetStart(x) ; i < GetEnd(x) ; ++i )
			if ( hash_table[i].y == y ) return i;
		return -1;
	}

	void initializeBIG( int n )
	{
		hash_table.re_allocate( HashDefaultRelativeSize * n );
		hash_table.m_num = hash_table.m_size;
		iVector<int> tmp;
		for ( int i = 0 ; i < n ; ++i )
		{
			tmp.push_back(i);
		}
		srand(time(0));
		for ( int i = 0 ; i < n ; ++i )
		{
			long long bigint = rand();
			bigint *= rand();
			bigint %= (n-i);
			hash_map.push_back( tmp[(int)bigint] );
			tmp[(int)bigint] = tmp[n-i-1];
		}
		m_num = n+1;
		m_start = new int[m_num];
		for ( int i = 0 ; i < m_num ; ++i )
		{
			m_start[i] = -1;
		}
		m_EdgeNum = 0;
		occupied = 0;
	}

	void free_mem()
	{
		delete[] m_start;
		hash_table.free_mem();
		hash_map.free_mem();
	}

	void initializeSMALL( int n )
	{
		hash_table.re_allocate( n );
		hash_table.m_num = n;
		m_start = new int[1];
		m_num = 1;
		m_EdgeNum = 0;
		occupied = 0;
	}

	void clean()
	{
		for ( int i = 0 ; i < hash_table.m_num ; ++i )
		{
			hash_table[i].x = CodeEdgeRange;
		}
		for ( int i = 0 ; i < m_num ; ++i )
		{
			m_start[i] = -1;
		}
		occupied = 0;
	}

	void SetAverageSpace( int r )
	{
		AverageSpace = hash_table.m_num/r;
	}

	int hash_functionBIG( int x , int y )
	{
		return HashDefaultRelativeSize*hash_map[x]+hash_map[y]%AverageSpace;
	}

	int hash_functionSMALL( int x , int y )
	{
		return y;
	}

	void Sort()
	{
		int i , j;
		for ( i = 0 , j = 0 ; i < hash_table.m_num ; ++i )
		{
			if ( hash_table[i].x < CodeEdgeRange )
			{
				if ( j != i )
				{
					hash_table[j] = hash_table[i];
					hash_table[i].x = CodeEdgeRange;
				}
				j++;
			}
		}
		hash_table.m_num = j;
		hash_table.Sort();
		hash_table.m_num = hash_table.m_size;
	}

	void insert( int x , int y , pWeight w , int code , bool bBIG )
	{
		int p;
		if ( bBIG )
		{
			p = hash_functionBIG( x , y );
		}
		else 
		{
			p = hash_functionSMALL( x , y );
			if ( occupied >= hash_table.m_num )
			{
				printf("FULL!!!!\n");
				printf("FULL!!!!\n");
				printf("FULL!!!!\n");
				printf("FULL!!!!\n");
				printf("!!!!!!!!\n");
			}
		}

		p %= hash_table.m_size;
		
		//find first x,y
		while ( hash_table[p].x != CodeEdgeRange && (hash_table[p].x != x || hash_table[p].y != y) )
		{
			p++;
			if ( p >= hash_table.m_size ) p -= hash_table.m_size;
		}

		if ( hash_table[p].x == CodeEdgeRange )
		{
			//first time x,y
			hash_table[p].x = x;
			hash_table[p].y = y;
			hash_table[p].weight = w;
			hash_table[p].code = code;
			hash_table[p].jump = 0;
			occupied++;
			return;
		}

		//update current x,y
		bool placed = false;
		int last_p = -1;
		int dis = 0;

		while ( code > 0 )
		{
			if ( w == hash_table[p].weight )
			{	
				if ( last_p != -1 )
				{
					hash_table[last_p].jump = dis;
				}
				hash_table[p].code |= code;
				placed = true;
				last_p = p;
				dis = 0;
			}
			else if ( hash_table[p].weight < w )
			{
				if ( last_p != -1 )
				{
					hash_table[last_p].jump = dis;
				}
				code &= (~hash_table[p].code);
				last_p = p;
				dis = 0;
			}
			else
			{
				hash_table[p].code &= (~code);
				if ( hash_table[p].code == 0 )
				{
					if ( !placed )
					{
						if ( last_p != -1 )
						{
							hash_table[last_p].jump = dis;
						}
						hash_table[p].weight = w;
						hash_table[p].code = code;						
						placed = true;
						last_p = p;
						dis = 0;
					}
					else
					{
						hash_table[p].x = CodeEdgeRange;
						occupied--;
					}
				}
				else
				{
					if ( last_p != -1 )
					{
						hash_table[last_p].jump = dis;
					}
					last_p = p;
					dis = 0;					
				}
			}
			
			if ( hash_table[p].jump == 0 ) break;
			else
			{
				dis += hash_table[p].jump;
				p += hash_table[p].jump;
				if ( p >= hash_table.m_size ) p -= hash_table.m_size;				
			}
		}

		//add new shortcut
		if ( !placed && code > 0 )
		{
			while ( hash_table[p].x != CodeEdgeRange )
			{
				p++;
				dis++;
				if ( p >= hash_table.m_size ) p -= hash_table.m_size;
			}
			
			hash_table[p].x = x;
			hash_table[p].y = y;
			hash_table[p].weight = w;
			hash_table[p].code = code;
			hash_table[p].jump = 0;
			hash_table[last_p].jump = dis;
			occupied++;
		}
		else
		{
			hash_table[last_p].jump = 0;
		}
	}

	void Indexing( bool bsorted )
	{
		//sort
		if ( !bsorted )
		{
			hash_table.Sort();
		}

		//indexing
		for ( int i = 0 ; i < hash_table.m_num ; ++i )
		{
			if ( i == 0 || hash_table[i].x != hash_table[i-1].x )
			{
				m_start[hash_table[i].x] = i;
			}
			if ( hash_table[i].x == CodeEdgeRange )
			{
				m_EdgeNum = i;
				break;
			}
		}

		for ( int i = m_num-1 ; i >= 0 ; --i )
		{
			if ( m_start[i] == -1 )
			{
				m_start[i] = m_start[i+1];
			}
		}
	}

	int GetStart( int x )
	{
		return m_start[x];
	}

	int GetEnd( int x )
	{
		return m_start[x+1];
	}
};

struct iHash
{
	int m_size;
	int* m_data;

	iHash( int n )
	{
		m_size = n;
		m_data = new int[m_size];
	}

	int find( int x )
	{
		int p = x % m_size;

		while ( m_data[p] > -1 )
		{
			if ( m_data[p] == x ) return p;
			p++;
			if ( p == m_size ) p = 0;
		}

		return p;
	}

	bool exist( int x )
	{
		return (m_data[find(x)] == x);
	}

	void insert( int x )
	{
		m_data[find(x)] = x;
	}

	void erase( int x )
	{
		m_data[find(x)] = -1;
	}

	void clean()
	{
		for ( int i = 0 ; i < m_size ; ++i )
			m_data[i] = -1;
	}
};

template <typename _T>
struct rBuffer
{
	_T* data;
	int ID, d_size, d_index, f_index, d_num, f_num;
	long long num;
	bool file_reverse, EoF;
	char file_name[32];
	FILE* file;
	
	rBuffer()
	{
	}

	void SmallBuffer( int id, int part, int size )
	{
		d_size = size;
		data = new _T[d_size];
		
		file_reverse = false;
		f_num = 1;
		f_index = 0;

		d_index = 0;		

		num = 0;

		EoF = false;

		sprintf(file_name, "%d.%d", id, part);
		file = fopen( file_name, "rb" );

		Refresh();
	}

	rBuffer( int id , int size, bool reverse, int file_num )
	{
		//each buffer must have a unique ID
		ID = id;
		//allocate
		d_size = size;
		data = new _T[d_size];

		file_reverse = reverse;
		f_num = file_num;

		if ( file_reverse ) f_index = f_num;
		else f_index = -1;

		d_index = d_num = 0;

		num = 0;

		EoF = false;

		file = NULL;
	}

	void Refresh()
	{
		d_num = fread( data, sizeof(_T), d_size, file );
		num += d_num;
	}

	bool NextFile()
	{
		if ( file_reverse ) f_index--;
		else f_index++;

		if ( f_index < 0 || f_index == f_num ) 
		{
			delete[] data;
			fclose(file);
			return false;
		}

		if ( file ) fclose(file);
		sprintf( file_name, "%d.%d", ID, f_index );
		file = fopen( file_name, "rb" );

		Refresh();

		return true;
	}

	bool GetNext( _T& x )
	{
		if ( d_index == d_num )
		{
			if ( d_num < d_size )
			{
				if ( !NextFile() ) 
				{
					EoF = true;
					return false;
				}
			}
			else
			{
				Refresh();
				if ( d_num == 0 )
				{
					if ( !NextFile() ) 
					{
						EoF = true;
						return false;
					}
				}
				else d_index = 0;
			}
		}

		x = data[d_index];
		d_index++;
		return true;
	}

	void GetCurrent( _T& x )
	{
		x = data[d_index];
		d_index++;

		if ( d_index == d_num )
		{
			if ( d_num < d_size )
			{
				if ( !NextFile() ) 
				{
					EoF = true;
					return;
				}
			}
			else
			{
				Refresh();
				if ( d_num == 0 )
				{
					if ( !NextFile() ) 
					{
						EoF = true;
						return;
					}
				}
				else d_index = 0;
			}
		}
	}

	_T top()
	{
		return data[d_index];
	}
};

template <typename _T>
struct wBuffer
{
	_T* data;
	int ID, d_size, d_index, f_num;
	long long num;
	bool b_unique, b_reverse, b_single;
	FILE* file;

	wBuffer( int id , int size, bool single, bool reverse, bool unique )
	{
		//each buffer must have a unique ID
		ID = id;
		//allocate
		d_size = size;
		data = new _T[d_size];

		b_single = single;
		b_reverse = reverse;
		b_unique = unique;

		f_num = 0;
		d_index = 0;
		num = 0;

		file = NULL;
	}

	void ToFile( bool EoF )
	{
		char file_name[32];
		sprintf( file_name, "%d.%d", ID, f_num );

		if ( b_single )
		{
			if ( !file ) file = fopen( file_name, "wb" );
		}
		else
		{
			file = fopen( file_name, "wb" );
			f_num++;
		}

		if ( b_reverse )
		{
			for ( int i = 0 ; i < d_index/2 ; ++i )
			{
				_T tmp = data[i];
				data[i] = data[d_index-1-i];
				data[d_index-1-i] = tmp;
			}
		}
		else if ( b_unique )
		{
			printf("sorting...");
			sort( data, data+d_index );
			printf("finished...");
			int j = 0;
			for ( int i = 0 ; i < d_index ; ++i )
				if ( !(data[i] == data[j]) )
				{
					++j;
					data[j] = data[i];
				}
			d_index = j+1;
			printf("unique finished\n");
		}

		fwrite( data, sizeof(_T), d_index, file );

		if ( d_size > 1024*1024 ) printf("wBuffer refresh: %d(out of %d) edges printed\n", d_index, d_size);

		num += d_index;

		if ( !b_single ) fclose(file);

		if ( EoF ) 
		{
			if ( b_single ) fclose(file);
			delete[] data;
		}
	}
	
	void push_back( _T x )
	{
		if ( d_index == d_size )
		{
			ToFile( false );
			d_index = 0;
		}
		data[d_index] = x;
		d_index++;
	}
};

struct oConfig
{
	FILE* file;

	oConfig()
	{
		file = fopen("config.txt","w");		
	}

	~oConfig()
	{
		fclose(file);
	}

	void AddConfig( string name, int setting )
	{
		fprintf( file, "%s %d\n", name.c_str(), setting );
	}
};

struct iConfig
{
	FILE* file;
	iVector<Key_Value<string, int> > configs;

	iConfig()
	{
		file = fopen("config.txt","r");

		char tmp[100];
		int setting;

		for ( ; fscanf( file, "%s %d", tmp, &setting ) == 2 ; )
		{
			configs.push_back( Key_Value<string,int>( string(tmp), setting ) );
		}

		fclose(file);
	}

	int LoadConfig( string name )
	{
		for ( int i = 0 ; i < configs.m_num ; ++i )
		{
			if ( configs[i].key == name ) return configs[i].value;
		}

		printf("No such config: %s\n", name.c_str());
		return -1;
	}
};	
struct TreeNode
{
	int left_cnt, cnt, left, right, father, depth;
};
struct LevelTree
{
	iVector<TreeNode> tree;//each node contains level and offset

	iVector<int> preorder;

	int N, node_num, depth_limit;

	int build_tree( int first, int last, int depth, int* l2n )
	{
		if ( first == last )
		{
			int p = l2n[first];
			tree[p].left = tree[p].right = -1;
			tree[p].left_cnt = tree[p].cnt = 1;
			tree[p].depth = depth;
			return p;
		}

		int mid = (first+last)/2;
		int p = node_num++;

		tree[p].depth = depth;

		int q = build_tree( first, mid, depth+1, l2n );
		tree[p].left = q;
		tree[q].father = p;
		
		tree[p].left_cnt = tree[q].cnt;

		q = build_tree( mid+1, last, depth+1, l2n );
		tree[p].right = q;
		tree[q].father = p;

		tree[p].cnt = tree[p].left_cnt + tree[q].cnt;

		return p;
	}

	void preorder_traverse( int x )
	{
		if ( tree[x].left != -1 ) preorder_traverse( tree[x].left );
		if ( x < N ) preorder.push_back( x );
		if ( tree[x].right != -1 ) preorder_traverse( tree[x].right );
	}

	void rebuild_tree()
	{
		printf("rebuilding tree...\n");
		preorder.clean();
		preorder_traverse( N );		
		node_num = N;
		build_tree( 0, N-1, 0, preorder.m_data );
	}

	void initialize( int* l2n, int n )
	{
		N = n;
		depth_limit = (int)(log((double)N)/log(2.0)*10);

		tree.re_allocate( 3*N );
		tree.m_num = 3*N;

		node_num = N;
		build_tree( 0, N-1, 0, l2n );
		tree[N].father = -1;

		preorder.clean();
		preorder_traverse( N );
		for ( int i = 0 ; i < N ; ++i )
		{
			printf("%d ",preorder[i]);
		}
		printf("\n");

		if ( node_num > 2*N ) printf("Building Tree Logic Error!!!\n");
	}

	int getLevel( int x )
	{
		int sum = 0;

		for ( int f ; x != N ; x = f )
		{
			f = tree[x].father;
			if ( tree[f].right == x ) sum += tree[f].left_cnt;
		}

		return sum;
	}

	void swap( int x , int y )
	{
		//delete x
		for ( int p = x, f ; p != N ; p = f )
		{
			f = tree[p].father;

			if ( tree[f].left == p )
			{
				tree[f].left_cnt--;

				if ( p == x )
				{
					tree[f].left = -1;
				}				
			}
			else if ( tree[f].right == p )
			{
				if ( p == x )
				{
					tree[f].right = -1;
				}
			}
			else
			{
				printf("Tree Structure Error, father: %d, child: %d\n",f,p);
			}

			tree[f].cnt--;
		}

		//insert x after y
		int f = tree[y].father;
		if ( tree[f].left == y )
		{
			tree[f].left = node_num++;
		}
		else if ( tree[f].right == y )
		{
			tree[f].right = node_num++;
		}
		else
		{
			printf("Tree Structure Error, father: %d, child: %d\n",f,y);
		}

		tree[node_num-1].left = y;
		tree[node_num-1].right = x;
		tree[node_num-1].father = f;
		tree[node_num-1].depth = tree[f].depth+1;
		tree[node_num-1].left_cnt = 1;
		tree[node_num-1].cnt = 2;

		tree[x].father = node_num-1;
		tree[x].depth = tree[node_num-1].depth+1;

		tree[y].father = node_num-1;
		tree[y].depth = tree[node_num-1].depth+1;

		if ( tree[x].depth > depth_limit )
		{
			rebuild_tree();
		}

		for ( int p = node_num-1 ; p != N ; p = f )
		{
			f = tree[p].father;
			if ( tree[f].left == p )
			{
				tree[f].left_cnt++;
			}			
			tree[f].cnt++;
		}
	}
};

struct LevelHash
{
	iVector<int> N2L, L2N;

	int C, R, Last;

	void initialize( int* level, int n , int c )
	{
		C = c;

		R = 0;

		Last = n*c+n;

		N2L.re_allocate( n );
		N2L.m_num = n;

		L2N.re_allocate( n*c+n+n );
		L2N.m_num = n*c+n+n;

		for ( int i = 0 ; i < n*c+n+n ; ++i )
		{
			L2N[i] = -1;
		}

		for ( int i = 0 ; i < n ; ++i )
		{
			N2L[i] = level[i]*c+n;

			L2N[level[i]*c+n] = i;
		}
	}

	void insert( int x, int y )
	{
		for ( int i = y-1 ; i > 0 ; --i )
		{
			if ( L2N[i] == -1 )
			{
				L2N[i] = x;
				N2L[x] = i;
				if ( y-i > 1000 )
				{
					reconstruct();
				}
				return;
			}
			else
			{
				int tmp = L2N[i];
				L2N[i] = x;
				N2L[x] = i;
				x = tmp;
			}
		}

		printf("BIG ERROR!!!!\n");
	}

	void reconstruct()
	{
		R++;
		for ( int i = 0 , j = 0; i < L2N.m_num ; ++i )
		{
			if ( L2N[i] != -1 )
			{
				N2L[L2N[i]] = N2L.m_num+j*C;
				j++;
				L2N[i] = -1;
			}
		}

		for ( int i = 0 ; i < N2L.m_num ; ++i )
		{
			if ( N2L[i] != -1 ) L2N[N2L[i]] = i;
		}

		Last = N2L.m_num*(C+1);
	}

	void remove( int x )
	{
		L2N[N2L[x]] = -1;
		N2L[x] = -1;
	}

	bool nottop( int k, int x )
	{
		return ( N2L[x] >= N2L.m_num+C*k );
	}

	void swap( int x , int y )
	{
		remove(x);

		if ( y == Last )
		{
			N2L[x] = y;
			L2N[y] = x;
			Last++;
		}
		else
		{		
			insert( x, y );
		}
	}
};
