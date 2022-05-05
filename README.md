# HGTphyloDetect
Horizontal gene transfer (HGT) refers to the exchange of genetic material
between disparate groups of organisms other than from parent to offspring,
which has been confirmed as a very significant factor in adaptive
evolution, disease emergence and metabolic shift that can act across various
species. However, current methods for HGT detection are usually not automatic,
narrow applicable or unavailable to use. In this work, we developed a versatile computational
toolbox named HGTphyloDetect by combining a high-throughput pipeline together with phylogenetic analysis
to facilitate comprehensive investigation of the potential mechanism for HGT events. Tests on two case
studies suggest that this approach could effectively
identify horizontally acquired genes with high accuracy. In-depth phylogenetic analysis further facilitates
the navigation of the potential donors and detailed gene transmission process. The HGTphyloDetect toolbox
is designed for ease of use and could accurately find HGT events with a very low false discovery rate
in a high-throughput manner.

## HGT identification pipeline
![image](https://github.com/SysBioChalmers/HGTphyloDetect/blob/master/doc/HGT_workflow.png)

## Example
We provide a user-friendly example for small test, users just need to prepare a Fasta file including protein id and protein sequence, note that protein id should be from the GenBank protein db. 
(1) If you are now in the HGTphyloDetect directory, just enter into the folder [example](https://github.com/SysBioChalmers/HGTphyloDetect/tree/master/example) via the command line:
```linux
cd example
```

(2) Users can run the script for the input file: 
```linux
python HGT_workflow.py input.fasta
```

(3) Then our software would generate the output results for this gene/protein, i.e., Alien index, E value and donor information:
|Gene/Protein | Alien index | E value | Donor id | 	Donor taxonomy |
| --------------- |:------------:|:------------:|:---------------:|:---------------:|
| AAT92670 |   199.17|    3.16e-87|  WP_208929673|  Bacteria/Firmicutes|

## Contact
* [Le Yuan](https://www.chalmers.se/en/Staff/Pages/leyu.aspx) ([@leyuan](https://github.com/le-yuan)), Chalmers University of Technology, Sweden
* [Eduard J. Kerkhoven](https://www.chalmers.se/en/staff/Pages/Eduard-Kerkhoven.aspx) ([@edkerk](https://github.com/edkerk)), Chalmers University of Technology, Sweden

