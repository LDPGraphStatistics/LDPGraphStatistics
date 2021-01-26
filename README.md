# LDPGraphStatistics

This is a source code of "Locally Differentially Private Analysis of Graph Statistics".

# Directory Structure
- cpp/			&emsp;C++ codes (put the required files under this directory; see cpp/README.md).
- data/			&emsp;Output data (obtained by running codes).
  - IMDB/			&emsp;Output data in the IMDB dataset.
  - Orkut/			&emsp;Output data in the Orkut dataset.
- python/		&emsp;Python codes.
- run_EvalLossLocal.sh	  &emsp;Shell script to evaluate the privacy and utility of graph algorithms.
- LICENSE.txt		&emsp;MIT license.
- README.md		&emsp;This file.

# Usage

**(1) Install**

Install StatsLib (see cpp/README.md).

Install C/C++ codes as follows.
```
$ cd cpp/
$ make
$ cd ../
```

**(2) Download and preprocess IMDB**

Download the [IMDB dataset](https://www.cise.ufl.edu/research/sparse/matrices/Pajek/IMDB.html) and place the dataset in data/IMDB/.


Run the following commands.

```
$ cd python/
$ python3 Read_IMDB.py data/IMDB/IMDB.mtx data/IMDB/edges.csv data/IMDB/deg.csv
$ cd ../
```

Then the edge file (edges.csv) and degree file (deg.csv) will be output in data/IMDB/.

**(3) Download and preprocess Orkut**

Download the [Orkut dataset](http://snap.stanford.edu/data/com-Orkut.html) and place the dataset in data/Orkut/.


Run the following commands.

```
$ cd python/
$ python3 Read_IMDB.py data/Orkut/com-orkut.ungraph.txt data/Orkut/edges.csv data/Orkut/deg.csv
$ cd ../
```

Then the edge file (edges.csv) and degree file (deg.csv) will be output in data/IMDB/.

**(4) Evaluate the privacy and utility of graph algorithms**

Run the following commands.

```
$ chmod +x run_EvalLossLocal.sh
$ ./run_EvalLossLocal.sh
```

Then experimental results for centralized DP (CentralLap_triangle, CentralLap_2star), one-round LDP (LocalRR_triangle, LocalLap_2star), and two-rounds LDP (Local2Rounds_triangle) graph algorithms (n=10000, epsilon=1) will be output in data/IMDB/ and data/Orkut/.

To change n and epsilon, change the 2nd and 3rd parameters in EvalLossLocal (for more details, see Usage of EvalLossLocal).

# Execution Environment
We used CentOS 7.5 with gcc 4.8.5 and python 3.6.5.

# External Libraries used in our source code.
- [StatsLib](https://www.kthohr.com/statslib.html) is distributed under the [Apache License 2.0](https://github.com/kthohr/stats/blob/master/LICENSE).
