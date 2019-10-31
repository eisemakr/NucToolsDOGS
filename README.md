# NucToolsDOGS
Additional scripts for the nuctools package by homeveg.
There are basically 2 scripts handling the output of compare2conditions.pl 
pos2reg.R converts the position-wise differences between 2 experimental conditions to .bed-file like regions. 
DOGS2Genes.R is a script making use of the Biomart package and shell promted bedtool intersect commands to gain insights into 
where the differentially occupied regions (DOGS = Differentially occupied genes) can be found.

Overlap_mouse_men.R finds orhtologous genes being differentially occupied in a mouse and a human sample. 
