# TIMES
This code accompanies the following papers:

* Inferring Temporal Information from a Snapshot of a Dynamic Network\
Jithin K. Sreedharan, Abram Magner, Ananth Grama, and Wojciech Szpankowski.\
_Nature Scientific Reports 2019_
* [TIMES: Temporal Information Maximally Extracted from Structures](https://dl.acm.org/citation.cfm?id=3186105)\
Abram Magner, Jithin K. Sreedharan, Ananth Grama, and Wojciech Szpankowski\
_ACM International Conference on World Wide Web (WWW) 2018_

#### Code
The code is written in C++. Use the `Makefile` and `Makefile.config` to make changes to the compilation process. It can be compiled either using g++7 or clang with support for C++11. The code make use of [SNAP library](https://snap.stanford.edu/snap/index.html) for graph data structures.

#### Data
Most of the data is taken from [SNAP database](https://snap.stanford.edu/data/index.html) with the exception of brain data which is taken from [Human Connectome Project](https://www.humanconnectome.org/study/hcp-young-adult/document/extensively-processed-fmri-data-documentation).
