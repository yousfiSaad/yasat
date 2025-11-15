# CI/CD Pipeline Documentation

## Overview

YASAT uses GitHub Actions for continuous integration, automated releases, and Docker image distribution. This document describes the CI/CD pipeline and how to use it.

## Table of Contents

1. [CI Pipeline](#ci-pipeline)
2. [Language Bindings](#language-bindings)
3. [Release Process](#release-process)
4. [Docker Images](#docker-images)
5. [Package Distribution](#package-distribution)
6. [Development Workflow](#development-workflow)

---

## CI Pipeline

### Workflows

The CI pipeline consists of three main workflows:

1. **CI Workflow** (`.github/workflows/ci.yml`) - Runs on every push and PR
2. **Release Workflow** (`.github/workflows/release.yml`) - Runs on version tags
3. **Docker Workflow** (`.github/workflows/docker.yml`) - Builds and publishes Docker images

### CI Workflow Details

**Triggers:**
- Push to `main`, `dev`, or `claude/**` branches
- Pull requests to `main` or `dev` branches

**Jobs:**

#### 1. Build and Test Matrix

Builds and tests on multiple platforms and configurations:

| Platform | Build Modes | Artifacts |
|----------|-------------|-----------|
| Ubuntu Latest | Release, Debug | Binary, Shared Library (.so), Static Library (.a) |
| macOS Latest | Release, Debug | Binary, Shared Library (.dylib), Static Library (.a) |
| Windows Latest | Release | Binary (.exe), DLL, Static Library (.a) |

**What it does:**
- Compiles YASAT for each platform
- Runs full test suite
- Builds shared and static libraries (release mode only)
- Uploads build artifacts for 7 days

#### 2. Code Quality Checks

- Runs `cppcheck` for static analysis
- Compiles with `-Werror` to catch all warnings

**What it tests (Linux/Release only):**
- C API tests: 24 tests
- Python binding tests: 25 tests
- Go binding tests: 17 tests
- Total: 82 tests across all language bindings

**View CI Results:**

Visit the Actions tab in your GitHub repository to see CI runs:
```
https://github.com/yousfiSaad/yasat/actions
```

---

## Language Bindings

YASAT provides a C API and bindings for multiple languages, all automatically tested in CI/CD.

### C API

Located in `src/c_api/`, provides a clean C interface to the SAT solver:

```c
#include <yasat.h>

// Create solver
yasat_solver* solver = yasat_solver_create();

// Add clauses (example: (1 OR 2) AND (NOT 1 OR 3))
yasat_add_clause(solver, (int[]){1, 2, 0}, 3);
yasat_add_clause(solver, (int[]){-1, 3, 0}, 3);

// Solve
yasat_result result = yasat_solve(solver);

// Clean up
yasat_solver_destroy(solver);
```

**C API Testing in CI:**
- 24 comprehensive tests
- Memory leak detection
- Edge case handling

### Python Bindings

Located in `bindings/python/`, provides a Pythonic interface:

```python
from yasat import Solver

solver = Solver()
solver.add_clause([1, 2])     # (1 OR 2)
solver.add_clause([-1, 3])    # (NOT 1 OR 3)

if solver.solve():
    print("SAT")
    print(solver.get_assignment())
else:
    print("UNSAT")
```

**Python Testing in CI:**
- 25 comprehensive tests
- Installation via `setup.py`
- Compatible with Python 3.7+

**Install from source:**
```bash
cd bindings/python
python setup.py install
```

### Go Bindings

Located in `bindings/go/`, provides idiomatic Go interface:

```go
import "github.com/yousfiSaad/yasat/bindings/go/yasat"

solver, _ := yasat.NewSolver()
defer solver.Close()

solver.AddClause([]int{1, 2})    // (1 OR 2)
solver.AddClause([]int{-1, 3})   // (NOT 1 OR 3)

if sat, _ := solver.Solve(); sat {
    fmt.Println("SAT")
    fmt.Println(solver.GetAssignment())
} else {
    fmt.Println("UNSAT")
}
```

**Go Testing in CI:**
- 17 comprehensive tests
- Full Go module support (`go.mod`)
- Compatible with Go 1.21+

**Install:**
```bash
go get github.com/yousfiSaad/yasat/bindings/go/yasat
```

### CI/CD Integration

All language bindings are:
- ✅ Automatically tested on every PR/push (Linux)
- ✅ Included in release packages
- ✅ Available in Docker images (`/usr/local/share/yasat/bindings/`)
- ✅ Documented with examples and tests

---

## Release Process

### Creating a Release

Releases are triggered automatically when you push a version tag:

```bash
# Tag a new version
git tag v1.0.0

# Push the tag to trigger release workflow
git push origin v1.0.0
```

### Version Tag Format

Use semantic versioning: `vMAJOR.MINOR.PATCH`

Examples:
- `v1.0.0` - Initial release
- `v1.1.0` - Minor update (new features, backward compatible)
- `v1.0.1` - Patch release (bug fixes)
- `v2.0.0` - Major update (breaking changes)

### What the Release Workflow Does

1. **Creates GitHub Release**
   - Extracts version from tag
   - Generates changelog from commits since last tag
   - Creates a new GitHub release

2. **Builds Release Artifacts** for all platforms:
   - Linux (x86_64): `yasat-linux-x86_64-VERSION.tar.gz`
   - macOS (x86_64): `yasat-macos-x86_64-VERSION.tar.gz`
   - Windows (x86_64): `yasat-windows-x86_64-VERSION.zip`

3. **Each package contains:**
   - Binary executable (`yasat` or `yasat.exe`)
   - Shared library (`.so`, `.dylib`, or `.dll`)
   - Static library (`.a`)
   - C++ header files (`include/` directory)
   - C API (`c_api/` directory)
   - Python bindings (`python/` directory)
   - Go bindings (`go/` directory)
   - README and license files
   - SHA256 checksum file

### Release Artifacts

All release artifacts are available on the GitHub Releases page:
```
https://github.com/yousfiSaad/yasat/releases
```

**Download and verify:**

```bash
# Download release
wget https://github.com/yousfiSaad/yasat/releases/download/v1.0.0/yasat-linux-x86_64-1.0.0.tar.gz

# Verify checksum
wget https://github.com/yousfiSaad/yasat/releases/download/v1.0.0/yasat-linux-x86_64-1.0.0.tar.gz.sha256
sha256sum -c yasat-linux-x86_64-1.0.0.tar.gz.sha256

# Extract
tar xzf yasat-linux-x86_64-1.0.0.tar.gz
cd yasat-linux-x86_64

# Use the binary
./yasat --help

# Or install the libraries
sudo cp libyasat.so* /usr/local/lib/
sudo cp -r include/* /usr/local/include/
sudo ldconfig
```

---

## Docker Images

### Available Images

Docker images are automatically built and published to GitHub Container Registry (ghcr.io).

**Image locations:**
- GitHub Container Registry: `ghcr.io/yousfisaad/yasat`
- Docker Hub: `yousfisaad/yasat` (requires configuration)

### Docker Image Tags

| Tag | Description | When Updated |
|-----|-------------|--------------|
| `latest` | Latest main branch build | Every push to `main` |
| `dev` | Latest dev branch build | Every push to `dev` |
| `vX.Y.Z` | Specific version | When version tag is pushed |
| `X.Y` | Major.minor version | When version tag is pushed |
| `X` | Major version | When version tag is pushed |
| `sha-XXXXXXX` | Specific commit | Every push |

### Using Docker Images

**Pull and run:**

```bash
# Pull latest image
docker pull ghcr.io/yousfisaad/yasat:latest

# Run YASAT
docker run --rm ghcr.io/yousfisaad/yasat:latest --help

# Solve a CNF file
docker run --rm -v $(pwd):/data ghcr.io/yousfisaad/yasat:latest /data/problem.cnf

# Use specific version
docker pull ghcr.io/yousfisaad/yasat:v1.0.0
docker run --rm ghcr.io/yousfisaad/yasat:v1.0.0 --version
```

**What's in the Docker image:**

- YASAT binary at `/usr/local/bin/yasat`
- Shared library at `/usr/local/lib/libyasat.so*`
- Static library at `/usr/local/lib/libyasat.a`
- Header files at `/usr/local/include/yasat/`
- Minimal Debian base (bookworm-slim)

**Using the shared library in Docker:**

```dockerfile
FROM ghcr.io/yousfisaad/yasat:latest as yasat

FROM debian:bookworm-slim

# Copy YASAT libraries from official image
COPY --from=yasat /usr/local/lib/libyasat.so* /usr/local/lib/
COPY --from=yasat /usr/local/include/yasat /usr/local/include/yasat

RUN ldconfig

# Your application code here
COPY myapp /usr/local/bin/myapp

CMD ["myapp"]
```

### Building Docker Images Locally

```bash
# Build image
docker build -t yasat:local .

# Run locally built image
docker run --rm yasat:local --help

# Test with a CNF file
docker run --rm -v $(pwd)/tests/cnf:/data yasat:local /data/simple_sat.cnf
```

---

## Package Distribution

### GitHub Releases (Primary)

**All platforms** are distributed via GitHub Releases:
- Pre-compiled binaries
- Shared and static libraries
- Header files
- SHA256 checksums

### Docker Hub (Optional)

To enable Docker Hub distribution:

1. Create Docker Hub account and repository
2. Add GitHub Secrets:
   - `DOCKERHUB_USERNAME`: Your Docker Hub username
   - `DOCKERHUB_TOKEN`: Docker Hub access token
3. Uncomment Docker Hub steps in `.github/workflows/docker.yml`

### Future Package Managers

#### Homebrew (macOS/Linux)

Create a Homebrew tap for easy installation:

```bash
brew tap yousfisaad/yasat
brew install yasat
```

**Setup:** Create a tap repository with a formula.

#### vcpkg (C++ Package Manager)

Add YASAT to vcpkg registry for C++ projects:

```bash
vcpkg install yasat
```

**Setup:** Create a vcpkg port with build instructions.

#### PyPI (Python Bindings)

Distribute Python bindings via PyPI:

```bash
pip install yasat
```

**Setup:** Create Python package with FFI bindings.

#### crates.io (Rust Bindings)

Publish Rust bindings to crates.io:

```toml
[dependencies]
yasat = "1.0"
```

**Setup:** Create Rust crate with FFI bindings.

---

## Development Workflow

### Branch Strategy

- `main` - Stable releases only
- `dev` - Development branch, merged to main for releases
- `claude/*` - Feature branches for automated development
- Feature branches - Use descriptive names

### CI in Pull Requests

When you open a PR:

1. **Automatic CI runs:**
   - Builds on all platforms
   - Runs full test suite
   - Runs code quality checks

2. **Required checks:**
   - All builds must pass
   - All tests must pass
   - Code quality checks must pass

3. **Review artifacts:**
   - Download build artifacts from the PR's CI run
   - Test binaries before merging

### Testing Locally Before Push

```bash
# Run full local build and test
make clean
make release
make debug
make lib
make test

# Check for warnings
make clean
CXXFLAGS="-std=c++17 -Wall -Wextra -Wpedantic -Werror" make release

# Test Docker build locally
docker build -t yasat:test .
docker run --rm yasat:test --help
```

### Release Checklist

Before creating a release tag:

- [ ] All CI checks pass on `main` branch
- [ ] Version number updated in `Makefile` (VERSION variable)
- [ ] Changelog reviewed (auto-generated from commits)
- [ ] All tests pass locally
- [ ] Docker image builds successfully
- [ ] Documentation is up to date

**Create release:**

```bash
# Update version in Makefile
vim Makefile  # Update VERSION = X.Y.Z

# Commit version bump
git add Makefile
git commit -m "Bump version to X.Y.Z"
git push origin main

# Create and push tag
git tag vX.Y.Z
git push origin vX.Y.Z

# Monitor release workflow
# Visit: https://github.com/yousfiSaad/yasat/actions
```

---

## Troubleshooting

### CI Failures

**Build fails on specific platform:**
- Check the CI logs for that platform
- Reproduce locally using the same compiler/OS
- Fix and push changes

**Tests fail:**
- Review test output in CI logs
- Run tests locally: `make test`
- Check for platform-specific issues

**Docker build fails:**
- Build locally: `docker build -t yasat:debug .`
- Check Dockerfile syntax
- Ensure all source files are included (check `.dockerignore`)

### Release Issues

**Release workflow fails:**
- Check that the tag follows `vX.Y.Z` format
- Ensure `GITHUB_TOKEN` has write permissions
- Review workflow logs for specific errors

**Missing artifacts:**
- Check if all platform builds completed successfully
- Retry the release workflow if needed
- Manually upload missing artifacts if necessary

### Docker Issues

**Image won't push:**
- Verify GitHub Container Registry permissions
- Check if `GITHUB_TOKEN` is available
- For Docker Hub, verify credentials are set

**Image too large:**
- Review `.dockerignore` to exclude unnecessary files
- Use multi-stage builds (already implemented)
- Check for unused dependencies

---

## CI/CD Configuration Files

### File Structure

```
.github/
└── workflows/
    ├── ci.yml          # Main CI pipeline
    ├── release.yml     # Release automation
    └── docker.yml      # Docker builds

Dockerfile              # Docker image definition
.dockerignore          # Docker build exclusions
Makefile               # Build system with library targets
```

### Key Makefile Targets

| Target | Description |
|--------|-------------|
| `make` | Build release binary (default) |
| `make release` | Build optimized binary |
| `make debug` | Build debug binary with symbols |
| `make lib` | Build both shared and static libraries |
| `make shared` | Build shared library only |
| `make static` | Build static library only |
| `make test` | Run test suite |
| `make clean` | Remove all build artifacts |
| `make install` | Install libraries and headers (requires sudo) |

---

## CI/CD Best Practices

### For Contributors

1. **Always run tests locally** before pushing
2. **Write clear commit messages** (they appear in changelog)
3. **Keep commits atomic** (one logical change per commit)
4. **Update tests** when adding features
5. **Document breaking changes** in commit messages

### For Maintainers

1. **Review CI results** before merging PRs
2. **Test release candidates** before tagging
3. **Keep workflows updated** with latest actions versions
4. **Monitor Docker image sizes**
5. **Update documentation** with CI/CD changes

---

## Support

For CI/CD issues:
1. Check workflow logs in GitHub Actions
2. Review this documentation
3. Open an issue with CI/CD logs attached
4. Tag issues with `ci/cd` label

**Useful links:**
- [GitHub Actions Documentation](https://docs.github.com/en/actions)
- [Docker Documentation](https://docs.docker.com/)
- [Semantic Versioning](https://semver.org/)

---

**Last Updated:** 2024
**Maintained by:** YOUSFI Saad
