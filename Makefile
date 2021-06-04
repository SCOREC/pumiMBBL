KOKKOS_PATH = /lore/vittav/Kokkos/kokkos
KOKKOS_DEVICES = "Cuda,OpenMP"
EXE_NAME = "pumi-test"

SRC_DIR = src
SRC = $(wildcard ${SRC_DIR}/*.cpp)

default: build
	echo "Start Build"

STATIC=pumiMBBLGPU.a

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
EXE = ${EXE_NAME}.cuda
KOKKOS_ARCH = "Turing75"
KOKKOS_CUDA_OPTIONS = "enable_lambda"
else
CXX = g++
EXE = ${EXE_NAME}.host
KOKKOS_ARCH = "BDW"
endif

CXXFLAGS = -O3
LINK = ${CXX}
LINKFLAGS =

DEPFLAGS = -M

OBJ = $(SRC:$(SRC_DIR)/%.cpp=%.o)
LIB =

include $(KOKKOS_PATH)/Makefile.kokkos

build: $(EXE)

$(STATIC): $(OBJ)
		@echo "[Link (Static)]"
		@ar rcs $@ $^

$(EXE): $(OBJ) $(STATIC) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(STATIC) $(KOKKOS_LIBS) $(LIB) -o $(EXE)
	rm *.o *.tmp *.h libkokkos.a

clean: kokkos-clean
	rm -f *.o *.host *.cuda *.a *.dat

# Compilation rules

%.o : $(SRC_DIR)/%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c -o $@ $<

test: $(EXE)
	./$(EXE)
