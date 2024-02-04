# YASAT - Yet Another SAT Solver

## Introduction

YASAT is a SAT solver that implements CDCL using only the C++ standard library and does not have any external dependency.


## Build

```
./build.sh
```

## Usage
```
cat ./data/php_2p_3h.cnf | ./build/yasat
```

### Output
```
SAT [1, 1, 0, 0, 0, 1]
```

