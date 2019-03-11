# anomaly_zone

This is just a python3 interface for the anomaly_finder.py script written
by CW Linkem (https://github.com/cwlinkem/anomaly_zone)

If you use this, please cite the original publication:
Linkem, C. W., Minin, V. N., & Leache, A. D. (2016). Detecting the anomaly zone 
in species trees and evidence for a misleading signal in higher-level skink phylogeny 
(Squamata: Scincidae). Systematic Biology, 65(3), 465-477.

## Installation
There is only one prerequisite which won't be available with your standard Python 3 installation:
* [DendroPy](https://dendropy.org/)

## Usage
View the help menu by calling anomaly_finder.py with the <-h> flag:

```
$ ./anomaly_finder.py -h

Usage:  ./anomaly_finder.py -t <treefile> -f <nexus or newick> 

Description: anomaly_finder.py by CW Linkem

	Arguments:
		-t		: Tree file 
		-f		: Format of tree file: "nexus" or "newick"
		-h		: Displays help menu
```
