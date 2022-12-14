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
![image](https://github.com/SysBioChalmers/HGTphyloDetect/blob/master/doc/HGT_pipeline.png)

## Installation
Install the latest version with:
```bash
$ git clone https://github.com/SysBioChalmers/HGTphyloDetect.git
$ cd HGTphyloDetect
$ pip install -r requirements.txt
```
Please note: this is sufficient to run the HGT detection functionality of HGTphyloDetect. Additional software dependencies and installation instructions are specified in the [User tutorial](https://github.com/SysBioChalmers/HGTphyloDetect/blob/master/User%20tutorial.pdf).

## Example
We provide a user-friendly example for small test, users just need to prepare a FASTA file including protein id and protein sequence, note that protein id should be from the GenBank protein database.<br><br>
(1) If you are now in the HGTphyloDetect directory, just enter into the folder [example](https://github.com/SysBioChalmers/HGTphyloDetect/tree/master/example) via the command line:
```linux
cd example
```

(2.1) Then users can run the script for the input file (default AI value = 45, out_pct = 0.90): 
```linux
python HGT_workflow.py input.fasta
```

(2.2) If users want to change the default values for the parameters used in the pipeline, e.g., AI value = 40, out_pct = 0.80, just reset the constant value and run the following: 
```linux
python HGT_workflow.py input.fasta AI=40 out_pct=0.80
```

(3) Finally, our software could generate the output results as a file under the folder [example](https://github.com/SysBioChalmers/HGTphyloDetect/tree/master/example) for this gene/protein. The output file includes some important information, i.e., Alien index, E value and donor information. For example:
|Gene/Protein | Alien index | E value | Donor id | 	Donor taxonomy |
|:---------------:|:---------------:|:---------------:|:---------------:|:---------------:|
| AAT92670 |   199.18|    3.15e-87|  WP_208929673|  Bacteria/Firmicutes|

## User Manual
Please check the documentation [User tutorial](https://github.com/SysBioChalmers/HGTphyloDetect/blob/master/User%20tutorial.pdf) for the user manual. Here you can download this file for your case studies.

## Contact
* [Le Yuan](https://www.chalmers.se/en/Staff/Pages/leyu.aspx) ([@leyuan](https://github.com/le-yuan)), Chalmers University of Technology, Sweden
* [Eduard J. Kerkhoven](https://www.chalmers.se/en/staff/Pages/Eduard-Kerkhoven.aspx) ([@edkerk](https://github.com/edkerk)), Chalmers University of Technology, Sweden

