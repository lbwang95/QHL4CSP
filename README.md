QHL: A Fast Algorithm for Exact Constrained Shortest Path Search on Road Networks
========================================================================

This repository contains the source codes and data for this paper. 

Usage
---------------

### Environment

g++ version: 8.4.0 

OS: Ubuntu/CentOS (Unix)

### Compilation

cd code

make

### Execution

In the folder code/,

./index [network name] [query path]

There is also a default setting where

network name=NY and query_path=q1

The "index" instruction will directly build all the index and run the 10 query sets, as stated in paper, stored in the folder "data/[network name]" . The network name could be NY, BAY, or COL. The output of the screen contains many axillary details, including the index time and memory costs. The results of the 10 query sets, including the query time, the number of hoplinks and path concatenations, will all be put in the file "Results" in the folder "data/[network name]".

### Data Description

The data used in the experiments are stored in data/. The folder contains three folders corresponding to the three networks.

In the folder of each network, it contains two files named "USA-road-t.[network name].gr" and "USA-road-d.[network name].gr" corresponding to the travel time and the distance of each edges. The descriptions can be found in the DIMACS challenge (http://www.diag.uniroma1.it//challenge9/download.shtml).

It also contains the 10 query set files named q1~q5 and r1~r5 corresponding to those stated in the paper. Each query set file contains 1000 lines. Each line contains three values of a CSP query which are s, t, and C.

The file named "random" contains 50000 random queries for building the QHL index. Each line is a random CSP query.

The file named "order.txt" is used to build the tree decomposition. One can generate it by using the function "genorder()" in the index.cpp.



