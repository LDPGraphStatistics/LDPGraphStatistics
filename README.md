# LDP Graph Statistics

This is a source code of the following paper:

Jacob Imola, Takao Murakami, Kamalika Chaudhuri, "Locally Differentially Private Analysis of Graph Statistics," Proceedings of the 30th USENIX Security Symposium (USENIX Security 2021), pp.983-1000, 2021.

Full paper: https://arxiv.org/abs/2010.08688

# Directory Structure
- cpp/			&emsp;C++ codes (put the required files under this directory; see cpp/README.md).
- data/			&emsp;Output data (obtained by running codes).
  - IMDB/			&emsp;Output data in the IMDB dataset.
  - Orkut/			&emsp;Output data in the Orkut dataset.
  - BAGraph-m10/			&emsp;Output data in the the BA (Barabasi-Albert) graph dataset (m=10).
  - BAGraph-m50/			&emsp;Output data in the the BA (Barabasi-Albert) graph dataset (m=50).
- python/		&emsp;Python codes.
- results/		&emsp;Experimental results.
  - IMDB/			&emsp;Experimental results in the IMDB dataset.
  - Orkut/			&emsp;Experimental results in the Orkut dataset.
  - BAGraph/		&emsp;Experimental results in the BA (Barabasi-Albert) graph datasets.
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

**(4) Generate the BA (Barabasi-Albert) graphs**

Run the following commands (m: attachment parameter).

```
$ cd python/
$ python3 BAGraph.py 1000000 10 data/BAGraph-m10/edges.csv data/BAGraph-m10/deg.csv
$ python3 BAGraph.py 1000000 50 data/BAGraph-m50/edges.csv data/BAGraph-m50/deg.csv
$ cd ../
```

Then the edge file (edges.csv) and degree file (deg.csv) will be output in data/BAGraph-m10/ and data/BAGraph-m50/.

**(5) Evaluate the privacy and utility of graph algorithms**

Run the following commands ([Dataset] is "IMDB", "Orkut", "BAGraph-m10", or "BAGraph-m50").

```
$ chmod +x run_EvalLossLocal.sh
$ ./run_EvalLossLocal.sh [Dataset (IMDB/Orkut/BAGraph-m10/BAGraph-m50)]
```

Then experimental results for centralized DP (CentralLap_triangle, CentralLap_2star), one-round LDP (LocalRR_triangle, LocalLap_2star), and two-rounds LDP (Local2Rounds_triangle) graph algorithms (n=10000, epsilon=1) will be output in data/[Dataset]/.

To change n and epsilon, change the 2nd and 3rd parameters in EvalLossLocal (for more details, see Usage of EvalLossLocal). Experimental results in our paper can be found in results/.

# Execution Environment
We used CentOS 7.5 with gcc 4.8.5 and python 3.6.5.

# External Libraries used in our source code.
- [StatsLib](https://www.kthohr.com/statslib.html) is distributed under the [Apache License 2.0](https://github.com/kthohr/stats/blob/master/LICENSE).
