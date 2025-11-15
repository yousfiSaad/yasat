# YASAT Docker Image
# Multi-stage build for minimal final image size

# Build stage
FROM gcc:13 AS builder

WORKDIR /build

# Copy source files
COPY src/ ./src/
COPY Makefile ./
COPY build.sh ./

# Create build directory
RUN mkdir -p build/objects build/lib_objects

# Build everything: binary, shared library, and static library
RUN make clean && \
    make release && \
    make lib

# Verify builds
RUN ls -lh build/ && \
    ldd build/yasat && \
    file build/yasat && \
    file build/libyasat.so* && \
    file build/libyasat.a

# Runtime stage - minimal image
FROM debian:bookworm-slim

LABEL maintainer="YOUSFI Saad <yousfi.saad@gmail.com>"
LABEL description="YASAT - Yet Another SAT Solver with shared library support"
LABEL version="1.0.0"

# Install runtime dependencies only
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libstdc++6 \
    && rm -rf /var/lib/apt/lists/*

# Create directories
RUN mkdir -p /usr/local/bin \
             /usr/local/lib \
             /usr/local/include/yasat

# Copy binary from builder
COPY --from=builder /build/build/yasat /usr/local/bin/yasat

# Copy libraries from builder
COPY --from=builder /build/build/libyasat.so* /usr/local/lib/
COPY --from=builder /build/build/libyasat.a /usr/local/lib/

# Copy headers from builder
COPY --from=builder /build/src/headers/*.h /usr/local/include/yasat/

# Create symlinks for shared library
RUN cd /usr/local/lib && \
    ln -sf libyasat.so.1.0.0 libyasat.so.1 && \
    ln -sf libyasat.so.1.0.0 libyasat.so && \
    ldconfig

# Add library path to LD_LIBRARY_PATH
ENV LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}

# Verify installation
RUN yasat --help || true && \
    ls -lh /usr/local/lib/libyasat* && \
    ls -lh /usr/local/include/yasat/

# Set working directory for user data
WORKDIR /data

# Default command shows help
ENTRYPOINT ["yasat"]
CMD ["--help"]
