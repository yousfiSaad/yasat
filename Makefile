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
LIB_OBJ_DIR = $(BUILD_DIR)/lib_objects
SRC_DIR = src
INSTALL_PREFIX ?= /usr/local

# Source files for CLI
CLI_SOURCES = $(SRC_DIR)/main.cpp $(SRC_DIR)/implementations/CDCL_solver.cpp
CLI_OBJECTS = $(OBJ_DIR)/main.o $(OBJ_DIR)/CDCL_solver.o
CLI_TARGET = $(BUILD_DIR)/yasat

# Source files for shared library
LIB_SOURCES = $(SRC_DIR)/implementations/CDCL_solver.cpp $(SRC_DIR)/c_api/yasat.cpp
LIB_OBJECTS = $(LIB_OBJ_DIR)/CDCL_solver.o $(LIB_OBJ_DIR)/yasat.o

# Detect platform for shared library extension
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    LIB_EXT = dylib
    LIB_FLAGS = -dynamiclib -install_name @rpath/libyasat.$(LIB_EXT)
else
    LIB_EXT = so
    LIB_FLAGS = -shared
endif

LIB_TARGET = $(BUILD_DIR)/libyasat.$(LIB_EXT)

# Default build mode (release)
MODE ?= release

# Set flags based on mode
ifeq ($(MODE),debug)
    CXXFLAGS += $(DEBUG_FLAGS)
else
    CXXFLAGS += $(RELEASE_FLAGS)
endif

# Phony targets
.PHONY: all clean debug release test help lib install uninstall

# Default target - build both CLI and library
all: release lib

# Release build (CLI only)
release:
	@$(MAKE) $(CLI_TARGET) MODE=release

# Debug build (CLI only)
debug:
	@$(MAKE) $(CLI_TARGET) MODE=debug

# Build shared library
lib:
	@$(MAKE) $(LIB_TARGET) MODE=release

# Build shared library in debug mode
lib-debug:
	@$(MAKE) $(LIB_TARGET) MODE=debug

# CLI target
$(CLI_TARGET): $(CLI_OBJECTS)
	@echo "Linking $(CLI_TARGET)..."
	$(CXX) $(LDFLAGS) -o $@ $^
	@echo "Build complete: $(CLI_TARGET)"

# Shared library target
$(LIB_TARGET): $(LIB_OBJECTS)
	@echo "Linking shared library $(LIB_TARGET)..."
	$(CXX) $(LIB_FLAGS) $(LDFLAGS) -o $@ $^
	@echo "Build complete: $(LIB_TARGET)"

# Compile CLI objects (without -fPIC)
$(OBJ_DIR)/main.o: $(SRC_DIR)/main.cpp
	@mkdir -p $(OBJ_DIR)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJ_DIR)/CDCL_solver.o: $(SRC_DIR)/implementations/CDCL_solver.cpp
	@mkdir -p $(OBJ_DIR)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Compile library objects (with -fPIC for position-independent code)
$(LIB_OBJ_DIR)/CDCL_solver.o: $(SRC_DIR)/implementations/CDCL_solver.cpp
	@mkdir -p $(LIB_OBJ_DIR)
	@echo "Compiling $< (PIC)..."
	$(CXX) $(CXXFLAGS) -fPIC -c -o $@ $<

$(LIB_OBJ_DIR)/yasat.o: $(SRC_DIR)/c_api/yasat.cpp
	@mkdir -p $(LIB_OBJ_DIR)
	@echo "Compiling $< (PIC)..."
	$(CXX) $(CXXFLAGS) -fPIC -c -o $@ $<

# Install targets
install: lib $(CLI_TARGET)
	@echo "Installing YASAT..."
	install -d $(INSTALL_PREFIX)/lib
	install -d $(INSTALL_PREFIX)/include
	install -d $(INSTALL_PREFIX)/bin
	install -m 755 $(LIB_TARGET) $(INSTALL_PREFIX)/lib/
	install -m 644 $(SRC_DIR)/c_api/yasat.h $(INSTALL_PREFIX)/include/
	install -m 755 $(CLI_TARGET) $(INSTALL_PREFIX)/bin/
	@echo "Installation complete."
	@echo "Library: $(INSTALL_PREFIX)/lib/libyasat.$(LIB_EXT)"
	@echo "Header:  $(INSTALL_PREFIX)/include/yasat.h"
	@echo "Binary:  $(INSTALL_PREFIX)/bin/yasat"

# Uninstall
uninstall:
	@echo "Uninstalling YASAT..."
	rm -f $(INSTALL_PREFIX)/lib/libyasat.so
	rm -f $(INSTALL_PREFIX)/lib/libyasat.dylib
	rm -f $(INSTALL_PREFIX)/include/yasat.h
	rm -f $(INSTALL_PREFIX)/bin/yasat
	@echo "Uninstall complete."

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -f $(CLI_TARGET) $(LIB_TARGET)
	rm -f $(OBJ_DIR)/*.o $(LIB_OBJ_DIR)/*.o
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
	@echo "  all         Build CLI and shared library (default)"
	@echo "  release     Build optimized CLI version"
	@echo "  debug       Build debug CLI version"
	@echo "  lib         Build shared library (release)"
	@echo "  lib-debug   Build shared library (debug)"
	@echo "  install     Install library, headers, and binary"
	@echo "  uninstall   Remove installed files"
	@echo "  clean       Remove build artifacts"
	@echo "  test        Build and run test suite"
	@echo "  help        Show this help message"
	@echo ""
	@echo "Variables:"
	@echo "  INSTALL_PREFIX  Installation prefix (default: /usr/local)"
	@echo ""
	@echo "Examples:"
	@echo "  make              # Build CLI and library"
	@echo "  make lib          # Build library only"
	@echo "  make debug        # Build debug CLI"
	@echo "  make install      # Install to /usr/local"
	@echo "  make install INSTALL_PREFIX=~/.local  # Install to home"
	@echo "  make clean        # Clean build artifacts"
	@echo "  make test         # Run tests"
