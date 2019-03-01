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
To build the algorithm of finding arrival order of a given graph in edge-list format:
```bash
make givengraph
```
To build the algorithm of finding arrival order of various random graph models:
```bash
make randomgraph
```

**Utilities for processing temporal graphs:**
```bash
make process_temporal_graph
```
An example run
```bash
./process_temporal_graph -i:"FB_wall_network.txt" -choice:1
```
This will process the graph file `FB_wall_network.txt` which has edges listed in `u v t` format per line (for node pair: `(u,v)`, the time of edge creation is `t`). It will output `predicted_rank.txt` with `u rank` format (arrival `rank` for node `u`) per line.

To build all the algorithms
```bash
make all
```
Running `make` first time takes time as it needs to compile the `SNAP` library.

#### Editing `Makefile` file for Mac systems
The default C++ compiler for Mac is `clang` and is invoked even if we use `g++` command. Instead, to use the GNU C++ compiler, install it via `brew install gcc` command, and change to `CC = g++-8` (8 is the version; replace it with the installed version) under `else ifeq ($(UNAME), Darwin)`.

## Data
Most of the data is taken from [SNAP database](https://snap.stanford.edu/data/index.html) with the exception of brain data which is taken from [Human Connectome Project](https://www.humanconnectome.org/study/hcp-young-adult/document/extensively-processed-fmri-data-documentation).
<!-- For the brain Connectome data, please download the brain connectome data from [here](link). The data is cleaned matrix version of the original human connectome project data. The code to clean the data is available here.
 -->
