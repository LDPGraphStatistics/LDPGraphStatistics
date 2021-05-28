#!/bin/bash -x

if [ $# -ne 1 ]; then
    echo "USAGE: run_EvalLossLocal.sh [Dataset]"
    exit 1
fi

cd cpp/
./EvalLossLocal ../data/${1}/edges.csv 10000 1 1 10 0
./EvalLossLocal ../data/${1}/edges.csv 10000 1 1 10 1
./EvalLossLocal ../data/${1}/edges.csv 10000 1 1 10 3
cd ../
