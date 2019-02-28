include Makefile.config

MAIN = vertex_ordering_general_directed_graph
# MAIN = vertex_ordering_random_graph
# MAIN = vertex_ordering_sequential_predict
# MAIN = vertex_ordering_copying_model
# MAIN = vertex_ordering_random_graph_exp
# MAIN = process_temporal_graph
# ————————————————

SRC_DIR = ./src
BUILD_DIR = ./build
# Snap-3.0 is required for finding P_uv forthe optmization; for compiling vertex_ordering_sequential_predict
SNAP_DIR = ./libraries/Snap-3.0
# SNAP_DIR = ./libraries/Snap-4.0

# DEPH = $(SRC_DIR)/vertex_ordering.hpp
DEPCPP = funs_vertex_ordering

all: $(MAIN)

# COMPILE

$(MAIN): $(SRC_DIR)/$(MAIN).cpp $(DEPH) $(SRC_DIR)/$(DEPCPP).cpp $(EXSNAP)/Snap.o
	$(CC) $(CXXFLAGS) -o $(BUILD_DIR)/$(MAIN) $(SRC_DIR)/$(MAIN).cpp $(SRC_DIR)/$(DEPCPP).cpp $(EXSNAP)/Snap.o -I$(EXSNAP) -I$(EXSNAPADV) -I$(EXGLIB) -I$(EXSNAPEXP) $(LDFLAGS) $(LIBS)

$(EXSNAP)/Snap.o:
	make -C $(EXSNAP)

clean:
	rm -f *.o  $(MAIN)  $(MAIN).exe
	rm -rf Debug Release
