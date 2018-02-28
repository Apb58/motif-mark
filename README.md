## Motif_Marker v1.1

by Adrian Bubie


**Motif Marker** is a python script used to identify short repeat/sequence motifs around intron-exon boundaries of given DNA sequences. Motifs are user defined, and translated into regular expression search terms using the IUPAC nomenclature standard. Motif locations returned by regex search are used to mark sub-sequence locations a 1-to-1 pixel to basepair graph of the intron-exon sequence. 

![Example Graph](exon_graphs.svg)

Note that Motif Marker recognizes exon sequences as using capitalized characters (ATCG) and intron sequences using lowercase (atcg), as per UCSC Genome Browser sequence download format.

Motif Marker can currently handle any number of intron-exon sequences, in FASTA format; however, note that only *one* intron-exon sequenceper FASTA sequence is graphed. If you would like to graph multiple exons in a gene, each exon must be split into its own FASTA entry (see fasta files in `/test` directory for examples).

### Downloads and Requirements

This program requires python v3.4+ to run.

To run the script, download the `Motif_Marker.py` executable. Navigating to download directoy, from the command line, enter:

```
./Motif_Marker.py -h
```

for instructions on how to run and pass in the required files.


(Current version: v1.1)
