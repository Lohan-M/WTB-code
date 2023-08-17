# WTB-code
This repository contains the code for a still unpublished paper titled "Reachability Queries on Dynamic Temporal Bipartite Graphs" and allows the reproduction of the experiments described in this paper.
This repository is still in progress, the code is only partial as the full code is still unpolished. We are working on it and on the documentation which is still minimal. We are sorry for any inconvenience caused by the lack of documentation.

## Datasets
The datasets used for the experiments can be directly downloaded on the website [Konect](http://konect.cc/networks/). For example, 'Wikiquote edits (it)' is one of the datasets used (as there are two of them on the website, check the number of edges to take the same as in the paper, here it would be the one with 703 139 edges).

The file needed is the 'out' file, (such as: 'out.edit-itwikiquote'). Delete the first lines (having a '%' at the start) from the file and delete the last line being empty. The dataset is then ready to use.

## Usage
After downloading the code and datasets, you can first run this command:
```c++
g++ -std=c++17 main.cpp new_core.cpp -o test -O3
```
Then, run this command:
```c++
./test ./path-to-dataset
```
With ./path-to-dataset the path to the dataset you want to run the experiments on.

The main function in the file new_core.cpp (at the end of the file) executes the test functions for the different operations. You can comment/uncomment the test functions to choose which operations you want to test. You can also change the value '1000' at line 1223 which controls the number of repetitions of SPRQs and SSRQs.
