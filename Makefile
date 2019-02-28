GLIB = glib-core
SNAP = snap-core
GLIBADV = glib-adv
SNAPADV = snap-adv
SNAPEXP = snap-exp
EXGLIB = $(SNAP_DIR)/$(GLIB)
EXSNAP = $(SNAP_DIR)/$(SNAP)
EXGLIBADV = $(SNAP_DIR)/$(GLIBADV)
EXSNAPADV = $(SNAP_DIR)/$(SNAPADV)
EXSNAPEXP = $(SNAP_DIR)/$(SNAPEXP)
UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
	CC = g++
	CXXFLAGS += -std=c++11 -Wall
	CXXFLAGS += -O3 -DNDEBUG -fopenmp
	# turn on for crash debugging, get symbols with <prog> 2>&1 | c++filt
	#CXXFLAGS += -g -rdynamic
	#CXXFLAGS += -ggdb
	# turn on for OpenMP
	#CXXOPENMP = -fopenmp
	LIBS += -lrt
else ifeq ($(UNAME), Darwin)
	CC = g++-8
	CXXFLAGS += -std=c++11 -Wall
	CXXFLAGS += -O3 -DNDEBUG
	# Enable the below flags and comment the above flag to use gdb
	# CXXFLAGS += -g3 -DNDEBUG
	CLANG := $(shell g++-7 -v 2>&1 | grep clang | cut -d " " -f 2)
	ifneq ($(CLANG), LLVM)
		CXXFLAGS += -fopenmp
	else
		CXXFLAGS += -DNOMP
	endif

else ifeq ($(shell uname -o), Cygwin)
	CC = g++
	CXXFLAGS += -Wall -D__STDC_LIMIT_MACROS
	CXXFLAGS += -O3 -DNDEBUG
	CXXOPENMP = -fopenmp
endif


MAIN = times_givengraph
SRC_DIR = ./src
BUILD_DIR = ./build
# Snap-3.0 is required for finding P_uv forthe optmization; for compiling vertex_ordering_sequential_predict
# SNAP_DIR = ./libraries/Snap-3.0
SNAP_DIR = ./libraries/Snap-4.0
DEPCPP = times_functions

all: givengraph randomgraph

givengraph: $(SRC_DIR)/times_givengraph.cpp $(SRC_DIR)/$(DEPCPP).cpp $(EXSNAP)/Snap.o | build
	$(CC) $(CXXFLAGS) -o $(BUILD_DIR)/times_givengraph $(SRC_DIR)/times_givengraph.cpp $(SRC_DIR)/$(DEPCPP).cpp $(EXSNAP)/Snap.o -I$(EXSNAP) -I$(EXSNAPADV) -I$(EXGLIB) -I$(EXSNAPEXP) $(LDFLAGS) $(LIBS)

randomgraph: $(SRC_DIR)/times_randomgraph.cpp $(SRC_DIR)/$(DEPCPP).cpp $(EXSNAP)/Snap.o | build
	$(CC) $(CXXFLAGS) -o $(BUILD_DIR)/times_randomgraph $(SRC_DIR)/times_randomgraph.cpp $(SRC_DIR)/$(DEPCPP).cpp $(EXSNAP)/Snap.o -I$(EXSNAP) -I$(EXSNAPADV) -I$(EXGLIB) -I$(EXSNAPEXP) $(LDFLAGS) $(LIBS)

build:
	mkdir -p $@

$(MAIN): $(SRC_DIR)/$(MAIN).cpp $(SRC_DIR)/$(DEPCPP).cpp $(EXSNAP)/Snap.o
	$(CC) $(CXXFLAGS) -o $(BUILD_DIR)/$(MAIN) $(SRC_DIR)/$(MAIN).cpp $(SRC_DIR)/$(DEPCPP).cpp $(EXSNAP)/Snap.o -I$(EXSNAP) -I$(EXSNAPADV) -I$(EXGLIB) -I$(EXSNAPEXP) $(LDFLAGS) $(LIBS)

$(EXSNAP)/Snap.o: | $(EXSNAP)/Snap.o
	make -C $(EXSNAP)

clean:
	rm -f *.o  $(MAIN)  $(MAIN).exe
	rm -rf Debug Release
