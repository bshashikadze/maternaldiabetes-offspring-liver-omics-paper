multi-omics analysis
================
BS
09/12/2022

Datasets required for this script can be found at:
1. data_1_log2 - is output of the [script - line 301](https://github.com/bshashikadze/maternaldiabetes-offspring-liver-omics-paper/blob/main/proteomics/proteomics.Rmd)
2. data_2_log2 - online supplementary tables, supplementary table 2 (sheet 1)
3. [Conditions file](https://github.com/bshashikadze/maternaldiabetes-offspring-liver-omics-paper/blob/main/clinical%20chemistry/Conditions.txt) (with information about samples and groups)

script of co-inertia analysis and visualization was adapted from
<https://github.com/bshashikadze/diabetes-lung-omics-paper/tree/main/multiomics%20coinertia>

CIA is based on: Meng, C., B. Kuster, A.C. Culhane, and A.M. Gholami, A
multivariate approach to the integration of multi-omics datasets. BMC
Bioinformatics, 2014. 15(1): p. 162.
<https://bioconductor.org/packages/release/bioc/html/omicade4.html>


