#!/bin/bash

for f in data/*.cnf
do
  echo -n "$f "
  cat "$f" | ./build/yasat
done

