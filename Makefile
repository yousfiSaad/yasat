# YASAT Makefile
# Build system for Yet Another SAT Solver

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -Wpedantic
LDFLAGS =

# Build modes
DEBUG_FLAGS = -g -DDEBUG -O0
RELEASE_FLAGS = -O3 -DNDEBUG

# Directories
BUILD_DIR = build
OBJ_DIR = $(BUILD_DIR)/objects
SRC_DIR = src

# Source files
SOURCES = $(SRC_DIR)/main.cpp $(SRC_DIR)/implementations/CDCL_solver.cpp
OBJECTS = $(OBJ_DIR)/main.o $(OBJ_DIR)/CDCL_solver.o
TARGET = $(BUILD_DIR)/yasat

# Default build mode (release)
MODE ?= release

# Set flags based on mode
ifeq ($(MODE),debug)
    CXXFLAGS += $(DEBUG_FLAGS)
else
    CXXFLAGS += $(RELEASE_FLAGS)
endif

# Phony targets
.PHONY: all clean debug release test help

# Default target
all: release

# Release build
release:
	@$(MAKE) $(TARGET) MODE=release

# Debug build
debug:
	@$(MAKE) $(TARGET) MODE=debug

# Main target
$(TARGET): $(OBJECTS)
	@echo "Linking $(TARGET)..."
	$(CXX) $(LDFLAGS) -o $@ $^
	@echo "Build complete: $(TARGET)"

# Compile main.o
$(OBJ_DIR)/main.o: $(SRC_DIR)/main.cpp
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Compile CDCL_solver.o
$(OBJ_DIR)/CDCL_solver.o: $(SRC_DIR)/implementations/CDCL_solver.cpp
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -f $(TARGET) $(OBJ_DIR)/*.o
	@echo "Clean complete."

# Run tests
test: release
	@echo "Running test suite..."
	@./tests/run_tests.sh

# Help target
help:
	@echo "YASAT Build System"
	@echo ""
	@echo "Targets:"
	@echo "  all         Build release version (default)"
	@echo "  release     Build optimized release version"
	@echo "  debug       Build debug version with symbols"
	@echo "  clean       Remove build artifacts"
	@echo "  test        Build and run test suite"
	@echo "  help        Show this help message"
	@echo ""
	@echo "Examples:"
	@echo "  make              # Build release version"
	@echo "  make debug        # Build debug version"
	@echo "  make clean        # Clean build artifacts"
	@echo "  make test         # Run tests"
