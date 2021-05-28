#!/bin/bash -x

cd cpp/
./EvalLossLocal ../data/IMDB/edges.csv 10000 1 1 10 0
./EvalLossLocal ../data/IMDB/edges.csv 10000 1 1 10 1
./EvalLossLocal ../data/IMDB/edges.csv 10000 1 1 10 3
./EvalLossLocal ../data/Orkut/edges.csv 10000 1 1 10 0
./EvalLossLocal ../data/Orkut/edges.csv 10000 1 1 10 1
./EvalLossLocal ../data/Orkut/edges.csv 10000 1 1 10 3
cd ../
