# TIMES
<!-- ![image](https://user-images.githubusercontent.com/19230005/53579065-e0fdf500-3b46-11e9-9437-6471a42f26f4.png) -->
<img src="https://user-images.githubusercontent.com/19230005/53579065-e0fdf500-3b46-11e9-9437-6471a42f26f4.png" width="600">


Given a single snapshot of a dynamic network, this code provides algorithms to infer the arrival order of the nodes.
It implements the optimal and approximate solutions of the problem. The optimal solution is a result of an integer programming formulation with coefficients found by random walk techniques.

The associated papers of this work are the following:
* [Inferring Temporal Information from a Snapshot of a Dynamic Network](https://rdcu.be/boQ5z)\
Jithin K. Sreedharan, Abram Magner, Ananth Grama, and Wojciech Szpankowski.\
_Nature Scientific Reports 2019_.\
Supplementary Material with all the mathematical details of the analysis and implementation is available [here](https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-019-38912-0/MediaObjects/41598_2019_38912_MOESM1_ESM.pdf)
* [TIMES: Temporal Information Maximally Extracted from Structures](https://dl.acm.org/citation.cfm?id=3186105)\
Abram Magner, Jithin K. Sreedharan, Ananth Grama, and Wojciech Szpankowski\
_ACM International Conference on World Wide Web (WWW) 2018_

## Getting Started
Download the repository or clone it with the following command
```bash
git clone https://github.com/jithin-k-sreedharan/times.git
```

### Prerequisites
- Compiler with C++11 support
- Gurobipy ([Gurobi](http://www.gurobi.com/) Python library). This is not required for the algorithms, but needed for running the script to find the optimal plot (Figure 3 in the papers)
- Python 2.7+ with Numpy and Scipy support.

### Installation
```bash
make
```
<!-- #### Code
The code is written in C++. Use the `Makefile` and `Makefile.config` to make changes to the compilation process. It can be compiled either using g++7 or clang with support for C++11. The code make use of [SNAP library](https://snap.stanford.edu/snap/index.html) for graph data structures.
 -->
## Data
Most of the data is taken from [SNAP database](https://snap.stanford.edu/data/index.html) with the exception of brain data which is taken from [Human Connectome Project](https://www.humanconnectome.org/study/hcp-young-adult/document/extensively-processed-fmri-data-documentation).
<!-- For the brain Connectome data, please download the brain connectome data from [here](link). The data is cleaned matrix version of the original human connectome project data. The code to clean the data is available here.
 -->
