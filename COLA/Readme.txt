Please cite our paper if you decide to user the code. The paper details are:

Sibo Wang, Xiaokui Xiao, Yin Yang, Wenqing Lin. 
Effective Indexing for Approximate Constrained Shortest Path Queries on Large Road Networks.
Proceedings of the VLDB Endowment (PVLDB), 10(2): 61-72, 2016.


In this project, we included the source code (COLA_code.zip) and the datasets, query sets we used in the experiments (COLA_datasets.zip). 

To complile
$cd COLA_zip/
$g++ cola.cpp  -O3 -o cola -std=c++11
You may use preprocessing.batch script to do the preprocessing

We have removed some datasets due to the size limitation of files in sourceforge.


To generate the partitions, we use the code provided by Yu Sun, which is publicly available at 
https://github.com/aldrichsun/Graph-Partitioning-with-Natural-Cuts


