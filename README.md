# TIMES
<!-- ![image](https://user-images.githubusercontent.com/19230005/53579065-e0fdf500-3b46-11e9-9437-6471a42f26f4.png) -->
<img src="https://user-images.githubusercontent.com/19230005/53579065-e0fdf500-3b46-11e9-9437-6471a42f26f4.png" width="600">


Given a single snapshot of a dynamic network, this code provides algorithms to infer the arrival order of the nodes.
It implements the optimal and approximate solutions of the problem. The optimal solution is a result of an integer programming formulation with coefficients found by random walk techniques.

The associated papers of this work are the following:
* [Inferring Temporal Information from a Snapshot of a Dynamic Network](https://rdcu.be/boQ5z)\
Jithin K. Sreedharan, Abram Magner, Ananth Grama, and Wojciech Szpankowski.\
_Nature Scientific Reports 2019_. [Supplementary Material](https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-019-38912-0/MediaObjects/41598_2019_38912_MOESM1_ESM.pdf) (details of analysis and implementation)
* [TIMES: Temporal Information Maximally Extracted from Structures](https://dl.acm.org/citation.cfm?id=3186105)\
Abram Magner, Jithin K. Sreedharan, Ananth Grama, and Wojciech Szpankowski\
_ACM International Conference on World Wide Web (WWW) 2018_

An example application from Scientific Reports paper:
<figure>
<img src="https://user-images.githubusercontent.com/19230005/53579668-fcb5cb00-3b47-11e9-8e39-dfd186865462.png" width="700">
<figcaption>
Human brain evolution deduced by our method from a network of human brain formed from Human Connectome Project (HCP) fMRI data.
</figcaption>
</figure>

## Getting Started
Download the repository or clone it with the following command
```bash
git clone https://github.com/jithin-k-sreedharan/times.git
```

### Prerequisites
- Compiler with C++11 support
- Gurobipy ([Gurobi](http://www.gurobi.com/) Python library). This is not required for the algorithms, but needed for running the script to find the optimal plot (Figure 3 in the papers)
- Python 2.7+ with Numpy and Scipy support.
- [SNAP](https://snap.stanford.edu/snap/index.html) graph library. The latest version of the library is already provided with the code and is available in `libraries` folder.

## Installation and Running the Algorithms
**Finding arrival order of a given graph in edge-list format:**
```bash
make times_givengraph
```
An example run:
```bash
# Find precision and recall of Peeling, Peeling+ and Maximum density precision-1 estimator
# Given Simple English Wikipedia graph edgelist and its node arrival order file
./times_givengraph -i:"../data/dynamic-simplewiki.txt" -ior:"../data/dynamic-simplewiki_data.csv" -choice:1
```
**Finding arrival order of various random graph models:**
```bash
make times_randomgraph
```
<!-- We use a generalized preferential attachment generator with parameters as follows
- `timen`: Total time-steps of the procedure
- `pralpha`: With this probability a new node will be added; with probbaility `(1-pralpha)`, new edges will be added between existing nodes.
- `vecp1`: Lower end of uniform distrbution for `m` (number of edges each new node brings into the graph) when a new node is added.
- `vecp2`: Upper end of uniform distrbution for `m` (number of edges each new node brings into the graph) when a new node is added.
- `prbeta`: With this probability end points of edges of the new node will be selected preferentially; with probability `(1-prbeta)` ebdpoints of edges of new node will be choosen uniformly at random.
- `vecq1`: Lower end of uniform distrbution for `m` (number of edges each new node brings into the graph) when edges between existing nodes are added;
- `vecq2`: Upper end of uniform distrbution for `m` (number of edges each new node brings into the graph) when edges between existing nodes are added.
- `prdelta`: With this probability,When adding edges between exisiting nodes, slource node be selected preferentially.
- `prgamma`: With this probability, when adding edges between exisiting nodes, Terminal node will be selected preferentially
 -->
 An example run:
```bash
# Find the average precision and recall of a generalized preferential attachment graph model
# See the beginning of /src/times_randomgraph.cpp for the explanation of parameters.
./times_randomgraph -timen:5000 -pralpha:0.75 -vecp1:5 -vecp2:50 -prbeta:0.5 -vecq1:5 -vecq2:50  -prdelta:0.5 -prgamma:0.5 -noruns:1000 -choice:0
```
**Utilities for processing temporal graphs:**
```bash
make process_temporal_graph
```
An example run:
```bash
# Process the graph file `FB_wall_network.txt` which has edges listed in `u v t` format per line
# i.e., for node pair: `(u,v)`, the time of edge creation is `t`.
# It will output `predicted_rank.txt` with `u rank` format (arrival `rank` for node `u`) per line.
./process_temporal_graph -i:"FB_wall_network.txt" -choice:1
```
**To reproduce the optimal inference plot**
1. Estimate p<sub>uv</sub>, the probability that for any two vertices u and v, vertex u is older than vertex v.
  ```bash
  make times_randomgraph_estimate_puv
  ```
  An example run:
  ```bash
  # See the beginning of /src/times_randomgraph_estimate_puv.cpp for the explanation of parameters.
  ./times_randomgraph_estimate_puv -timen:50 -pralpha:0.5 -vecp1:1 -vecp2:8 -prbeta:0.5 -vecq1:1 -vecq2:8  -prdelta:0.5 -prgamma:0.5 -noruns:100 -norunsMC:100 -choice:3
  ```
2. Run the IPython Jupyter file `src/times_optimization.ipynb`

**To build all the algorithms**
```bash
make all
```
Running `make` first time takes time as it needs to compile the `SNAP` library.

#### Compiling issues and solutions
- `Makefile` in Mac systems: The default C++ compiler for Mac is `clang` and is invoked even if we use `g++` command. Instead, to use the GNU C++ compiler, install it via `brew install gcc` command, and change to `CC = g++-8` (8 is the version; replace it with the installed version) under `else ifeq ($(UNAME), Darwin)`.
- `Snap` library issues an error when running `times_randomgraph_estimate_puv`: Line 723 of `Snap-4.0/snap-core/network.cpp` needs to be commented out to remove the error. It is already changed if you use the `Snap` library provided in the `libraries` folder.

## Data
Most of the data are taken from [SNAP database](https://snap.stanford.edu/data/index.html) with the exception of brain network data. We describe the procedure to collect raw brain data and to form brain networks out of it in a separate [repository](https://github.com/jithin-k-sreedharan/data_human_brain_networks).

## Acknowledgments
This work was supported by NSF Center for Science of Information (CSoI) Grant CCF-0939370, NSF Grants CCF-1524312 and CSR-1422338, and by NIH Grant 1U01CA198941-01.
<!-- For the brain Connectome data, please download the brain connectome data from [here](link). The data is cleaned matrix version of the original human connectome project data. The code to clean the data is available here.
 -->
