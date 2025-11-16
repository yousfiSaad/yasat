# FFI Callback Challenges: A Comprehensive Guide

## Table of Contents

1. [Introduction](#introduction)
2. [Why Callbacks Are Difficult Across FFI](#why-callbacks-are-difficult-across-ffi)
3. [Platform-Specific Challenges](#platform-specific-challenges)
4. [Solutions and Design Patterns](#solutions-and-design-patterns)
5. [Threading and Synchronization](#threading-and-synchronization)
6. [Complete Implementation Examples](#complete-implementation-examples)
7. [Recommendations for YASAT](#recommendations-for-yasat)
8. [References](#references)

---

## Introduction

When building a shared library (like YASAT) that will be called from multiple programming languages, one common requirement is supporting **callbacks** - functions written in the host language (Go, Python, Rust, etc.) that are called by the C/C++ library during execution.

For example, you might want:
- A progress callback that reports solving progress: `onProgress(iteration, conflicts)`
- A solution callback that receives satisfying assignments: `onSolution(variables)`
- A conflict callback for debugging: `onConflict(clause, level)`

**The key question**: Is implementing callbacks across FFI boundaries straightforward?

**The answer**: **No, callbacks are NOT straightforward** and come with significant complexity and platform-specific challenges.

This document explains why callbacks are difficult and provides practical solutions.

---

## Why Callbacks Are Difficult Across FFI

### 1. Calling Convention Mismatches

Different languages use different **calling conventions** (how arguments are passed, who cleans up the stack, register usage):

- **C/C++**: Uses system ABI (cdecl, stdcall on Windows)
- **Go**: Uses its own custom calling convention (Go ABI)
- **Python**: Uses C calling convention via ctypes
- **Rust**: Can use C calling convention with `extern "C"`

**Problem**: A C library expects a function pointer with C calling convention. If you pass a Go function pointer directly, the calling convention mismatch causes crashes.

### 2. Runtime Differences

Different language runtimes have fundamentally different execution models:

- **C/C++**: Direct machine code, no runtime overhead
- **Go**: Garbage collected, cooperative scheduling, goroutines
- **Python**: Interpreted (or JIT), garbage collected, Global Interpreter Lock (GIL)
- **Rust**: No runtime, but needs to respect ownership/lifetime rules

**Problem**: When C code calls back into Go/Python, it's crossing runtime boundaries. The host runtime may need to be "initialized" or "locked" for safe execution.

### 3. Memory Management and Lifetimes

Who owns the callback function? How long does it live?

- **C/C++**: Function pointers are just addresses, no lifetime tracking
- **Go**: Functions might be closures capturing variables; GC can move data
- **Python**: Functions are objects that can be garbage collected
- **Rust**: Callbacks must respect Rust's lifetime and borrow checker rules

**Problem**: The C library stores a function pointer, but the host language might deallocate or move the underlying function/closure.

### 4. Thread Safety

SAT solving is often single-threaded, but what if callbacks spawn work?

- **Python**: GIL must be held/released correctly
- **Go**: Callback might be called from a C thread, not a goroutine
- **Rust**: Must ensure thread-safety with `Send`/`Sync` traits

**Problem**: Calling back into a different language's runtime from an arbitrary thread can cause deadlocks or crashes.

---

## Platform-Specific Challenges

### Go: The Most Challenging Case

Go has unique challenges that make callbacks particularly difficult:

#### 1. **No Stable Function Addresses**

In Go, you **cannot** pass a Go function directly as a C function pointer:

```go
// THIS DOES NOT WORK
func myCallback(x int) {
    fmt.Println("Progress:", x)
}

// ❌ Cannot convert Go function to C function pointer
cCallback := C.progress_callback(myCallback)  // COMPILE ERROR
```

**Why?**: Go functions don't have stable addresses in memory. The Go compiler/runtime may move code around.

#### 2. **Closures and Captured Variables**

Go callbacks are often closures that capture variables:

```go
counter := 0
callback := func(x int) {
    counter += x  // Captures 'counter'
    fmt.Println("Total:", counter)
}
```

**Problem**: The captured variables are Go-managed memory. C has no way to access or preserve this context.

#### 3. **Goroutine Runtime Dependency**

Go functions assume they're running in a goroutine with access to the Go runtime:

- Allocations use Go's garbage collector
- Defer statements, panic/recover mechanisms
- Channel operations, select statements

**Problem**: If C calls the callback from a C thread (not initialized by Go), the Go runtime isn't available, causing crashes.

#### 4. **The Solution: Gateway Functions**

Go provides `//export` to create C-callable functions:

```go
//export GoProgressCallback
func GoProgressCallback(handle C.int, progress C.int) {
    // This is callable from C
}
```

**But**: These exported functions:
- Cannot be closures (no captured variables)
- Must have C-compatible signatures
- Require manual context management via handles

---

### Python: Moderate Difficulty

Python is easier than Go but still has challenges:

#### 1. **GIL (Global Interpreter Lock)**

Python's GIL must be held when executing Python code:

```c
// In C library calling Python callback
PyGILState_STATE gstate = PyGILState_Ensure();
callback(progress);  // Now safe to call Python
PyGILState_Release(gstate);
```

**Problem**: If you don't acquire the GIL, Python code execution will crash or corrupt memory.

#### 2. **Reference Counting**

Python uses reference counting for memory management:

```python
def make_callback():
    local_data = [1, 2, 3]
    def callback(x):
        print(sum(local_data) + x)
    return callback
```

**Problem**: The callback object and its captured `local_data` must be kept alive while C holds the function pointer.

#### 3. **ctypes.CFUNCTYPE**

Python provides `ctypes.CFUNCTYPE` to create C-callable function pointers:

```python
import ctypes

# Define callback signature
CALLBACK_TYPE = ctypes.CFUNCTYPE(None, ctypes.c_int)

# Create callback
def my_callback(x):
    print(f"Progress: {x}")

c_callback = CALLBACK_TYPE(my_callback)

# Pass to C library
lib.yasat_set_progress_callback(solver, c_callback)

# IMPORTANT: Keep reference to c_callback!
# If it gets garbage collected, C will have a dangling pointer
```

---

### Rust: Relatively Straightforward

Rust is the easiest case for callbacks:

```rust
extern "C" fn progress_callback(progress: i32) {
    println!("Progress: {}", progress);
}

unsafe {
    yasat_set_progress_callback(solver, Some(progress_callback));
}
```

**Why it works**:
- Rust functions with `extern "C"` use C calling convention
- No runtime to worry about
- Explicit control over lifetimes and thread safety

**Challenge**: Closures still require manual context management (similar to Go).

---

## Solutions and Design Patterns

### Pattern 1: Simple Same-Thread Callbacks (Recommended for YASAT)

**Design**: Keep callbacks simple, synchronous, and same-thread.

```c
// C API
typedef void (*yasat_progress_callback)(int iteration, int conflicts);

void yasat_set_progress_callback(yasat_solver* solver,
                                  yasat_progress_callback callback,
                                  void* user_data);
```

**Properties**:
- Callback is called from the same thread that called `yasat_solve()`
- Callback executes synchronously (blocking)
- No threading complexity
- User provides `user_data` for context

**Go Usage**:

```go
// Global registry to map handles to Go callbacks
var callbackRegistry = make(map[int]func(int, int))
var callbackCounter = 0
var callbackMutex sync.Mutex

//export GoProgressCallback
func GoProgressCallback(handle C.int, iteration C.int, conflicts C.int) {
    callbackMutex.Lock()
    callback := callbackRegistry[int(handle)]
    callbackMutex.Unlock()

    if callback != nil {
        callback(int(iteration), int(conflicts))
    }
}

// User-facing API
func (s *Solver) SetProgressCallback(callback func(iteration, conflicts int)) {
    callbackMutex.Lock()
    handle := callbackCounter
    callbackCounter++
    callbackRegistry[handle] = callback
    callbackMutex.Unlock()

    C.yasat_set_progress_callback(
        s.handle,
        C.yasat_progress_callback(C.GoProgressCallback),
        unsafe.Pointer(uintptr(handle)),
    )
}
```

**Python Usage**:

```python
import ctypes

PROGRESS_CALLBACK = ctypes.CFUNCTYPE(None, ctypes.c_void_p, ctypes.c_int, ctypes.c_int)

class Solver:
    def __init__(self):
        self.lib = ctypes.CDLL("libyasat.so")
        self.handle = self.lib.yasat_solver_create()
        self._callback = None  # Keep reference!

    def set_progress_callback(self, callback):
        def wrapper(user_data, iteration, conflicts):
            callback(iteration, conflicts)

        self._callback = PROGRESS_CALLBACK(wrapper)
        self.lib.yasat_set_progress_callback(
            self.handle,
            self._callback,
            None
        )

# Usage
def on_progress(iteration, conflicts):
    print(f"Iteration {iteration}: {conflicts} conflicts")

solver = Solver()
solver.set_progress_callback(on_progress)
solver.solve()
```

### Pattern 2: Global Registry Pattern (for Go Closures)

**Problem**: Go doesn't allow closures in callbacks.

**Solution**: Use a global map to associate integer handles with Go closures.

```go
var progressCallbacks = make(map[int]func(int, int))
var nextHandle = 1
var mutex sync.Mutex

//export ProgressCallbackGateway
func ProgressCallbackGateway(handle C.int, iteration C.int, conflicts C.int) {
    mutex.Lock()
    callback := progressCallbacks[int(handle)]
    mutex.Unlock()

    if callback != nil {
        callback(int(iteration), int(conflicts))
    }
}

func (s *Solver) SetProgressCallback(cb func(int, int)) {
    mutex.Lock()
    handle := nextHandle
    nextHandle++
    progressCallbacks[handle] = cb
    mutex.Unlock()

    C.yasat_set_progress_callback(
        s.cSolver,
        C.yasat_progress_callback(C.ProgressCallbackGateway),
        unsafe.Pointer(uintptr(handle)),
    )
}
```

**Cleanup** (important):

```go
func (s *Solver) Close() {
    if s.callbackHandle != 0 {
        mutex.Lock()
        delete(progressCallbacks, s.callbackHandle)
        mutex.Unlock()
    }
    C.yasat_solver_destroy(s.cSolver)
}
```

### Pattern 3: Polling Alternative (No Callbacks)

If callbacks are too complex, offer a polling-based API instead:

```c
// C API
typedef struct {
    int iteration;
    int conflicts;
    int decisions;
    bool is_solving;
} yasat_status;

yasat_status yasat_get_status(yasat_solver* solver);
```

**Go Usage**:

```go
go func() {
    ticker := time.NewTicker(100 * time.Millisecond)
    defer ticker.Stop()

    for range ticker.C {
        status := solver.GetStatus()
        if !status.IsSolving {
            break
        }
        fmt.Printf("Progress: %d iterations, %d conflicts\n",
                   status.Iteration, status.Conflicts)
    }
}()

result := solver.Solve()  // Blocks
```

**Advantages**:
- No callback complexity
- Each language can poll in its native way
- No cross-runtime issues

**Disadvantages**:
- Polling overhead
- Not real-time
- Requires thread-safe status updates

---

## Threading and Synchronization

### Challenge: Cross-Thread Callbacks

If your C library is multithreaded and calls callbacks from worker threads:

#### Python GIL Complications

```c
// In C worker thread
void worker_thread() {
    // ... solving work ...

    // MUST acquire GIL before calling Python
    PyGILState_STATE gstate = PyGILState_Ensure();
    progress_callback(iteration, conflicts);
    PyGILState_Release(gstate);
}
```

**Problem**: Acquiring the GIL from arbitrary threads is complex and can cause deadlocks if not done carefully.

#### Go Runtime Complications

```go
//export ProgressCallback
func ProgressCallback(iteration C.int) {
    // This might be called from a C thread!
    // Go runtime might not be initialized on this thread

    // Solution: Use runtime.LockOSThread() or ensure
    // C library calls us from the same thread as yasat_solve()
}
```

### Recommendation: Same-Thread Callbacks

For YASAT (currently single-threaded), **always call callbacks from the same thread** that called `yasat_solve()`. This avoids all threading complexity:

```c
bool yasat_solve(yasat_solver* solver) {
    // ... solving loop ...

    if (solver->progress_callback && iteration % 1000 == 0) {
        // Called from same thread as yasat_solve()
        solver->progress_callback(solver->user_data, iteration, conflicts);
    }

    // ... continue solving ...
}
```

---

## Complete Implementation Examples

### C API Design

```c
// yasat.h
#ifndef YASAT_H
#define YASAT_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct yasat_solver yasat_solver;

// Callback types
typedef void (*yasat_progress_callback)(void* user_data,
                                        int iteration,
                                        int conflicts);

typedef void (*yasat_solution_callback)(void* user_data,
                                        int num_vars,
                                        const int* assignment);

// Solver management
yasat_solver* yasat_solver_create(void);
void yasat_solver_destroy(yasat_solver* solver);

// Input
int yasat_add_clause(yasat_solver* solver, int* literals, int num_literals);
int yasat_parse_cnf_file(yasat_solver* solver, const char* filename);

// Callbacks
void yasat_set_progress_callback(yasat_solver* solver,
                                 yasat_progress_callback callback,
                                 void* user_data);

void yasat_set_solution_callback(yasat_solver* solver,
                                 yasat_solution_callback callback,
                                 void* user_data);

// Solving
typedef enum {
    YASAT_RESULT_SAT = 10,
    YASAT_RESULT_UNSAT = 20,
    YASAT_RESULT_UNKNOWN = 0,
    YASAT_RESULT_ERROR = -1
} yasat_result;

yasat_result yasat_solve(yasat_solver* solver);

// Solution retrieval (alternative to callback)
int yasat_get_assignment(yasat_solver* solver, int var);

#ifdef __cplusplus
}
#endif

#endif // YASAT_H
```

### C++ Implementation

```cpp
// yasat_c_api.cpp
#include "yasat.h"
#include "CDCL_solver.h"

struct yasat_solver {
    CDCL_solver cpp_solver;
    yasat_progress_callback progress_cb;
    void* progress_user_data;
    yasat_solution_callback solution_cb;
    void* solution_user_data;
};

yasat_solver* yasat_solver_create(void) {
    try {
        return new yasat_solver{};
    } catch (...) {
        return nullptr;
    }
}

void yasat_solver_destroy(yasat_solver* solver) {
    delete solver;
}

void yasat_set_progress_callback(yasat_solver* solver,
                                 yasat_progress_callback callback,
                                 void* user_data) {
    if (solver) {
        solver->progress_cb = callback;
        solver->progress_user_data = user_data;
    }
}

yasat_result yasat_solve(yasat_solver* solver) {
    if (!solver) return YASAT_RESULT_ERROR;

    try {
        // Modify CDCL_solver to call progress callback
        // (requires adding callback support to C++ code)
        bool result = solver->cpp_solver.solve();
        return result ? YASAT_RESULT_SAT : YASAT_RESULT_UNSAT;
    } catch (...) {
        return YASAT_RESULT_ERROR;
    }
}
```

### Go Bindings with Global Registry

```go
package yasat

/*
#cgo LDFLAGS: -L. -lyasat
#include "yasat.h"

extern void GoProgressCallback(void* handle, int iteration, int conflicts);
*/
import "C"
import (
    "sync"
    "unsafe"
)

// Global callback registry
var (
    progressCallbacks = make(map[int]func(int, int))
    nextCallbackHandle = 1
    callbackMutex sync.Mutex
)

//export GoProgressCallback
func GoProgressCallback(handle unsafe.Pointer, iteration C.int, conflicts C.int) {
    handleInt := int(uintptr(handle))

    callbackMutex.Lock()
    callback := progressCallbacks[handleInt]
    callbackMutex.Unlock()

    if callback != nil {
        callback(int(iteration), int(conflicts))
    }
}

type Solver struct {
    cSolver        *C.yasat_solver
    callbackHandle int
}

func NewSolver() (*Solver, error) {
    cSolver := C.yasat_solver_create()
    if cSolver == nil {
        return nil, errors.New("failed to create solver")
    }
    return &Solver{cSolver: cSolver}, nil
}

func (s *Solver) SetProgressCallback(callback func(iteration, conflicts int)) {
    callbackMutex.Lock()
    handle := nextCallbackHandle
    nextCallbackHandle++
    progressCallbacks[handle] = callback
    callbackMutex.Unlock()

    s.callbackHandle = handle

    C.yasat_set_progress_callback(
        s.cSolver,
        C.yasat_progress_callback(C.GoProgressCallback),
        unsafe.Pointer(uintptr(handle)),
    )
}

func (s *Solver) Solve() (bool, error) {
    result := C.yasat_solve(s.cSolver)

    switch result {
    case C.YASAT_RESULT_SAT:
        return true, nil
    case C.YASAT_RESULT_UNSAT:
        return false, nil
    default:
        return false, errors.New("solver error")
    }
}

func (s *Solver) Close() error {
    // Clean up callback
    if s.callbackHandle != 0 {
        callbackMutex.Lock()
        delete(progressCallbacks, s.callbackHandle)
        callbackMutex.Unlock()
    }

    C.yasat_solver_destroy(s.cSolver)
    return nil
}
```

### Python Bindings with ctypes

```python
import ctypes
from typing import Callable, Optional
from enum import IntEnum

# Load library
lib = ctypes.CDLL("libyasat.so")

# Define types
class YasatResult(IntEnum):
    SAT = 10
    UNSAT = 20
    UNKNOWN = 0
    ERROR = -1

# Define callback types
PROGRESS_CALLBACK = ctypes.CFUNCTYPE(
    None,                    # return type
    ctypes.c_void_p,        # user_data
    ctypes.c_int,           # iteration
    ctypes.c_int            # conflicts
)

# Define function signatures
lib.yasat_solver_create.restype = ctypes.c_void_p
lib.yasat_solver_destroy.argtypes = [ctypes.c_void_p]
lib.yasat_set_progress_callback.argtypes = [
    ctypes.c_void_p,        # solver
    PROGRESS_CALLBACK,       # callback
    ctypes.c_void_p         # user_data
]
lib.yasat_solve.argtypes = [ctypes.c_void_p]
lib.yasat_solve.restype = ctypes.c_int

class Solver:
    def __init__(self):
        self.handle = lib.yasat_solver_create()
        if not self.handle:
            raise RuntimeError("Failed to create solver")

        # CRITICAL: Keep reference to callback to prevent GC
        self._progress_callback = None

    def set_progress_callback(self, callback: Callable[[int, int], None]):
        """Set progress callback. Callback receives (iteration, conflicts)."""

        # Wrapper to handle user_data parameter
        def wrapper(user_data, iteration, conflicts):
            callback(iteration, conflicts)

        # Create C callback (and keep reference!)
        self._progress_callback = PROGRESS_CALLBACK(wrapper)

        # Register with C library
        lib.yasat_set_progress_callback(
            self.handle,
            self._progress_callback,
            None  # user_data
        )

    def solve(self) -> Optional[bool]:
        """Solve the CNF formula. Returns True (SAT), False (UNSAT), or None (error)."""
        result = lib.yasat_solve(self.handle)

        if result == YasatResult.SAT:
            return True
        elif result == YasatResult.UNSAT:
            return False
        else:
            return None

    def __del__(self):
        if hasattr(self, 'handle') and self.handle:
            lib.yasat_solver_destroy(self.handle)

# Usage example
def main():
    def on_progress(iteration, conflicts):
        if iteration % 1000 == 0:
            print(f"Iteration {iteration}: {conflicts} conflicts")

    solver = Solver()
    solver.set_progress_callback(on_progress)

    # ... add clauses ...

    result = solver.solve()
    if result is True:
        print("SAT")
    elif result is False:
        print("UNSAT")
    else:
        print("ERROR")
```

---

## Recommendations for YASAT

### 1. Start Simple: Progress Callback Only

For the initial shared library implementation, provide **one simple callback**:

```c
typedef void (*yasat_progress_callback)(void* user_data,
                                        int iteration,
                                        int conflicts);
```

**Properties**:
- Called every N iterations (e.g., every 1000)
- Same-thread execution (from within `yasat_solve()`)
- User provides context via `user_data` pointer

### 2. Use Global Registry Pattern for Go

Since Go is a primary target language, implement the global registry pattern shown above. Document this pattern clearly for users.

### 3. Provide Polling Alternative

Also provide a polling-based API for users who don't want callback complexity:

```c
typedef struct {
    int iteration;
    int conflicts;
    int decisions;
    int restarts;
} yasat_stats;

yasat_stats yasat_get_current_stats(yasat_solver* solver);
```

### 4. Document Reference Lifetime Requirements

Clearly document in Python bindings that users **must** keep references to callback objects:

```python
class Solver:
    def set_progress_callback(self, callback):
        # Store reference to prevent garbage collection
        self._callback_ref = CALLBACK_TYPE(callback)
        lib.yasat_set_callback(self.handle, self._callback_ref)
```

### 5. Keep Callbacks Synchronous

Don't support async/concurrent callbacks initially. Keep them simple and blocking:

- Callback executes on same thread
- Callback blocks solver until it returns
- No threading/GIL/goroutine complexity

### 6. Consider Timeout Alternative

Instead of interruption callbacks, provide timeout parameter:

```c
yasat_result yasat_solve_with_timeout(yasat_solver* solver,
                                      double timeout_seconds);
```

This avoids needing an interruption callback mechanism.

---

## References

### Real-World Examples

**SQLite** - C library with simple callback design:
```c
int sqlite3_exec(
    sqlite3*,                    /* Database handle */
    const char *sql,             /* SQL statement */
    int (*callback)(void*,int,char**,char**),  /* Callback */
    void *user_data,             /* User data for callback */
    char **errmsg                /* Error message */
);
```

**libgit2** - Progress callbacks with user data:
```c
typedef struct {
    int (*sideband_progress)(const char *str, int len, void *payload);
    int (*completion)(unsigned int completion_type, void *payload);
    void *payload;
} git_remote_callbacks;
```

### Documentation Resources

- **CGo Documentation**: https://pkg.go.dev/cmd/cgo
- **Python ctypes**: https://docs.python.org/3/library/ctypes.html
- **Rust FFI Guide**: https://doc.rust-lang.org/nomicon/ffi.html
- **Go Wiki on cgo**: https://github.com/golang/go/wiki/cgo

### Key Takeaways

1. **Callbacks across FFI are NOT straightforward** - they require careful design
2. **Go is the hardest** - requires global registry pattern for closures
3. **Python is moderate** - requires GIL awareness and reference management
4. **Rust is easiest** - direct `extern "C"` support
5. **Start simple** - single progress callback, same-thread, synchronous
6. **Provide alternatives** - polling APIs for users who want to avoid callbacks
7. **Document clearly** - especially lifetime/reference requirements

---

**For YASAT Development**:

When implementing the shared library, start with:
1. One simple progress callback
2. Clear examples for Go (with registry pattern)
3. Clear examples for Python (with reference management)
4. Thorough documentation of requirements and limitations

This will provide 80% of the value with 20% of the complexity.
