# Compiler and flags for Linux
CXX = g++
CXXFLAGS =  -std=c++1y -O3 -ffast-math -march=native
CXXFLAGS += -fopenmp
# If MKL is installed, comment first line and uncomment MKL directives
#LDFLAGS += -llapack -lblas
LDFLAGS += -Wl,--no-as-needed -liomp5 -lpthread -lm -ldl -llapack -larpack
### MKL directives
MKLFLAGS += -L${MKLROOT}/lib/intel64
MKLFLAGS += -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core 

# Directories
SRC_DIR = ../src
BUILD_DIR = .
INCLUDE_DIR = $(SRC_DIR)

# Source files and object files
SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRC))

# Header files (for dependency tracking)
HEADERS = $(wildcard $(SRC_DIR)/*.hpp)

# Target executable
TARGET = ../main

# Default target
all: $(TARGET)

# Link the executable
$(TARGET): $(OBJ)
	$(CXX) $(LDFLAGS) $(MKLFLAGS) -o $@ $^

# Compile object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -c $< -o $@

# Clean up build files
clean:
	rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/$(TARGET)
