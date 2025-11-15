#!/bin/bash

# YASAT Build Script
# Simple wrapper around the Makefile for easy building

# Create build directories if they don't exist
mkdir -p build/objects

# Build using the Makefile (defaults to release mode)
make "$@"
