# File and R codes used to investigate enrichment of top cis-sQTLs and top cis-eQTLs

## File downloaded from ENCODE (under “ENCODE_file” folder):
  * "eCLIP_url_208.txt"
  * "experiment_report_208.tsv"
  * "experiment_report_208_encode3.tsv"
  * "experiment_report_208_encode4.tsv"



## CODE:
_**DIVIDE sQTLs and eQTLs into 20 GROUPS**_ 

Group criteria with index:
| 	          | MAF <= 0.05	| (0.05, 0.1]	| (0.1, 0.2]  | (0.2, 0.3]  | (0.3, 0.5] |
|-----------------|-------------|---------------|-------------|-------------|------------|
| dis <= 1e+3	  | 1	        |  5	        |  9	      |  13	    |   17	 |
| (1e+3, 1e+4]	  | 2	        |  6	        |  10	      |  14	    |   18	 |
| (1e+4, 1e+5]	  | 3	        |  7	        |  11	      |  15	    |   19	 |
| (1e+5, 1e+6]	  | 4	        |  8	        |  12	      |  16	    |   20	 |

**s1_group_eqtl_sqtl.R**
  * Aim: group top sQTLs and top eQTLs
  * Output: "top_sqtl_eqtl_group.RData"

**s2_group_null_eqtl.R**
  * Aim: group null cis-eQTLs from eQTL analysis
  * Output: "null_eqtl_group.RData"

**s3_group_null_sqtl.R**
  * Aim: group null cis-sQTLs from sQTL analysis
  * Output: "null_sqtl_group.RData"


_**ANNOTATE TOP eQTLs and sQTLs to RBPs**_ 

**s4_address_exp_report.R**
  * Aim: address information of 208 eCLIP experiments from ENCODE: peak file, cell line, RBP, ENCODE version, url
  * Output: "dat_experiment_208.csv"

**s5_annot_qtl_eclip.R**
  * Aim: annotate eQTLs or sQTLs with RNA binding proteins (RBPs) from ENCODE.
  * Output: "match_eqtl_all_eclip208_sum.RData"; "match_sqtl_all_eclip208_sum.RData"

**s6_summary_annot_qtl.R**
  * Aim: calculated number of top sQTLs or eQTLs annotated to the corresponding RBP
  * Output: "annot_eqtl_all_freq.csv"; "annot_sqtl_all_freq.csv"


_**ANNOTATE NULL eQTLs and sQTLs to RBPs**_

**s7_annot_sqtlp4.R**
  * Aim: annotate distinct cis-sQTLs with RBPs from ENCODE
  * Output: "nullsqtl_annot_i.RData"

**s8_annot_eqtlp4.R**
  * Aim: annotate distinct cis-eQTLs with RBPs from ENCODE
  * Output: "nulleqtl_annot_i.RData"

**s9_summary_annot_nullqtl.R**
  * Aim: In previous steps, we annotate null sQTLs and null eQTLs after dividing them to 10000 SNPs per file. Here, we merge annotation outputs to two files (i.e., one for null sQTLs and one for eQTLs).
  * Output: "annot_sqtlp4.RData"; "annot_eqtl_p4.RData"


_**RANDOM SAMPLING**_

**s10_random_sampling_1000.R**
  * Aim: random sampling 1000 times for top sQTLs or top eQTLs based on the number of top QTLs in each of 20 groups
  * Output: "sampling.RData"

**s11_summary_sampling_annot.R**
  * Aim: annotate SNPs in each random sampling; calculate the number of null SNPs annotated for each RBP in each random sampling; calculate the average proportion of null SNPs annotated for each RBP across random sampling.
  * Output: "eqtl_rbp_sampling.csv"; "sqtl_rbp_sampling.csv"


_**ENRICHMENT VIOLIN PLOT**_

**s12_enrich_violin_plot.R**
  * Aim: calculate fold for sQTL and eQTL; print violin plot
  * Output: "plot_eclip.pdf"

