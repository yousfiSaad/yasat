# YASAT Makefile
# Build system for Yet Another SAT Solver

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -Wpedantic
LDFLAGS =

# Library version
VERSION = 1.0.0
VERSION_MAJOR = 1

# Build modes
DEBUG_FLAGS = -g -DDEBUG -O0
RELEASE_FLAGS = -O3 -DNDEBUG

# Directories
BUILD_DIR = build
OBJ_DIR = $(BUILD_DIR)/objects
LIB_OBJ_DIR = $(BUILD_DIR)/lib_objects
SRC_DIR = src
INCLUDE_DIR = $(SRC_DIR)/headers

# Source files
SOLVER_SRC = $(SRC_DIR)/implementations/CDCL_solver.cpp
MAIN_SRC = $(SRC_DIR)/main.cpp

# Object files
MAIN_OBJ = $(OBJ_DIR)/main.o
SOLVER_OBJ = $(OBJ_DIR)/CDCL_solver.o
SOLVER_PIC_OBJ = $(LIB_OBJ_DIR)/CDCL_solver.o

# Binary target
TARGET = $(BUILD_DIR)/yasat

# Library targets
LIB_NAME = yasat
STATIC_LIB = $(BUILD_DIR)/lib$(LIB_NAME).a

# Platform-specific shared library names
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    SHARED_LIB = $(BUILD_DIR)/lib$(LIB_NAME).so.$(VERSION)
    SHARED_LIB_LINK = $(BUILD_DIR)/lib$(LIB_NAME).so
    SHARED_LIB_SONAME = lib$(LIB_NAME).so.$(VERSION_MAJOR)
    SHARED_FLAGS = -shared -Wl,-soname,$(SHARED_LIB_SONAME)
endif
ifeq ($(UNAME_S),Darwin)
    SHARED_LIB = $(BUILD_DIR)/lib$(LIB_NAME).$(VERSION).dylib
    SHARED_LIB_LINK = $(BUILD_DIR)/lib$(LIB_NAME).dylib
    SHARED_FLAGS = -dynamiclib -install_name @rpath/lib$(LIB_NAME).$(VERSION_MAJOR).dylib
endif
ifneq (,$(findstring MINGW,$(UNAME_S)))
    SHARED_LIB = $(BUILD_DIR)/$(LIB_NAME).dll
    SHARED_LIB_LINK = $(SHARED_LIB)
    SHARED_FLAGS = -shared
endif
ifneq (,$(findstring MSYS,$(UNAME_S)))
    SHARED_LIB = $(BUILD_DIR)/$(LIB_NAME).dll
    SHARED_LIB_LINK = $(SHARED_LIB)
    SHARED_FLAGS = -shared
endif

# Default build mode (release)
MODE ?= release

# Set flags based on mode
ifeq ($(MODE),debug)
    CXXFLAGS += $(DEBUG_FLAGS)
else
    CXXFLAGS += $(RELEASE_FLAGS)
endif

# Phony targets
.PHONY: all clean debug release test help lib shared static install

# Default target
all: release

# Release build
release:
	@$(MAKE) $(TARGET) MODE=release

# Debug build
debug:
	@$(MAKE) $(TARGET) MODE=debug

# Build all libraries
lib: shared static

# Build shared library
shared:
	@$(MAKE) $(SHARED_LIB) MODE=release

# Build static library
static:
	@$(MAKE) $(STATIC_LIB) MODE=release

# Build directories
$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)

$(LIB_OBJ_DIR):
	@mkdir -p $(LIB_OBJ_DIR)

# Main binary target
$(TARGET): $(OBJ_DIR) $(MAIN_OBJ) $(SOLVER_OBJ)
	@echo "Linking $(TARGET)..."
	$(CXX) $(LDFLAGS) -o $@ $(MAIN_OBJ) $(SOLVER_OBJ)
	@echo "Build complete: $(TARGET)"

# Shared library target
$(SHARED_LIB): $(LIB_OBJ_DIR) $(SOLVER_PIC_OBJ)
	@echo "Building shared library $(SHARED_LIB)..."
	$(CXX) $(SHARED_FLAGS) -o $@ $(SOLVER_PIC_OBJ)
ifneq ($(SHARED_LIB),$(SHARED_LIB_LINK))
	@cd $(BUILD_DIR) && ln -sf $(notdir $(SHARED_LIB)) $(notdir $(SHARED_LIB_LINK))
endif
	@echo "Shared library built: $(SHARED_LIB)"

# Static library target
$(STATIC_LIB): $(LIB_OBJ_DIR) $(SOLVER_PIC_OBJ)
	@echo "Building static library $(STATIC_LIB)..."
	ar rcs $@ $(SOLVER_PIC_OBJ)
	@echo "Static library built: $(STATIC_LIB)"

# Compile main.o (no PIC needed for executable)
$(MAIN_OBJ): $(MAIN_SRC) | $(OBJ_DIR)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Compile CDCL_solver.o (no PIC for executable)
$(SOLVER_OBJ): $(SOLVER_SRC) | $(OBJ_DIR)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Compile CDCL_solver.o with PIC (for shared library)
$(SOLVER_PIC_OBJ): $(SOLVER_SRC) | $(LIB_OBJ_DIR)
	@echo "Compiling $< (PIC)..."
	$(CXX) $(CXXFLAGS) -fPIC -c -o $@ $<

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -f $(TARGET) $(OBJ_DIR)/*.o $(LIB_OBJ_DIR)/*.o
	rm -f $(SHARED_LIB) $(SHARED_LIB_LINK) $(STATIC_LIB)
	rm -f $(BUILD_DIR)/*.so* $(BUILD_DIR)/*.dylib* $(BUILD_DIR)/*.dll $(BUILD_DIR)/*.a
	@echo "Clean complete."

# Run tests
test: release
	@echo "Running test suite..."
	@./tests/run_tests.sh

# Install libraries and headers
PREFIX ?= /usr/local
INSTALL_LIB_DIR = $(PREFIX)/lib
INSTALL_INCLUDE_DIR = $(PREFIX)/include/yasat

install: lib
	@echo "Installing YASAT libraries and headers..."
	install -d $(INSTALL_LIB_DIR)
	install -d $(INSTALL_INCLUDE_DIR)
	install -m 644 $(SHARED_LIB) $(INSTALL_LIB_DIR)
	install -m 644 $(STATIC_LIB) $(INSTALL_LIB_DIR)
ifneq ($(SHARED_LIB),$(SHARED_LIB_LINK))
	cd $(INSTALL_LIB_DIR) && ln -sf $(notdir $(SHARED_LIB)) $(notdir $(SHARED_LIB_LINK))
endif
	install -m 644 $(INCLUDE_DIR)/*.h $(INSTALL_INCLUDE_DIR)
	@echo "Installation complete."
	@echo "Libraries installed to: $(INSTALL_LIB_DIR)"
	@echo "Headers installed to: $(INSTALL_INCLUDE_DIR)"

# Uninstall
uninstall:
	@echo "Uninstalling YASAT..."
	rm -f $(INSTALL_LIB_DIR)/lib$(LIB_NAME).*
	rm -rf $(INSTALL_INCLUDE_DIR)
	@echo "Uninstall complete."

# Help target
help:
	@echo "YASAT Build System"
	@echo ""
	@echo "Targets:"
	@echo "  all         Build release binary (default)"
	@echo "  release     Build optimized release binary"
	@echo "  debug       Build debug binary with symbols"
	@echo "  lib         Build both shared and static libraries"
	@echo "  shared      Build shared library (.so/.dylib/.dll)"
	@echo "  static      Build static library (.a)"
	@echo "  clean       Remove all build artifacts"
	@echo "  test        Build and run test suite"
	@echo "  install     Install libraries and headers (requires sudo)"
	@echo "  uninstall   Uninstall libraries and headers (requires sudo)"
	@echo "  help        Show this help message"
	@echo ""
	@echo "Examples:"
	@echo "  make              # Build release binary"
	@echo "  make debug        # Build debug binary"
	@echo "  make lib          # Build shared and static libraries"
	@echo "  make shared       # Build only shared library"
	@echo "  make clean        # Clean build artifacts"
	@echo "  make test         # Run tests"
	@echo "  sudo make install # Install to /usr/local"
