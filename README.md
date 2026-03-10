🔍 Overview of Analyses  
  
1. Download genome files from NCBI  
   These SLURM scripts automate downloading, verifying checksums, and decompressing genome FASTA files on an HPC cluster.  
  ・01_wget_genomes_parallel_protostome.slurm  
  ・02_compare_checksums_parallel_protostome.slurm  
  ・03_decompress_fna.slurm (with 03_decompress_fna.sh)  
  
2. Calculate dinucleotide counts, frequencies, Observed/Expected ratios for whole-genomes  
   This pipeline calculates dinucleotide counts for full genome assemblies using sliding window analyses.  
  ・04_whole_genome_calc  
  
3. Calculate dinucleotide counts, frequencies, Observed/Expected ratios for BUSCO genes using Miniprot as a gene prediction tool  
   Gene prediction and extraction of BUSCO regions using Miniprot, followed by motif frequency calculations  
     
   Primary pipeline (recommended):  
  ・05_BUSCOs_with_Miniprot_pipeline  
  
   Alternative pipelines:  
    ・06_BUSCOs_with_Augustus_pipeline  
    ・07_all_genes_with_Augustus_pipeline  
       These two pipelines are almost the same as 05_BUSCOs_with_Miniprot_pipeline but using different gene prediction tools.  
       06_BUSCOs_with_Augustus_pipeline is the pipeline using BUSCO genes with Augustus, and  
       07_all_genes_with_Augustus_pipeline is the pipeline for extracting all genes in a genome with Augustus.  
       For studies analyzing lifespan-associated nucleotide motif evolution, we recommend using Miniprot BUSCO genes  
       as this is the fastest method and they can capture sufficient lifespan information as the other two gene sets.    
   
4. Removed  directory
   ・08_######
   This directory was removed during code revision and is no longer included.
  
5. Generate supplementary figures and tables in the published article (additional analyses)  
   These scripts reproduce main figures and supplementary plots/tables in the manuscript  
  ・09_upstream_sliding_window_figures  
  ・10_downstream_sliding_window_figures  
  ・11_BUSCO_Miniprot_BUSCO_Augustus_allgenes_comparison_figures  
  ・Some figures and tables were generated from code in 05_10_dinucleotide_calc_pipeline_excluding_5spp directory
  
  
💻 Script Types and Execution  
SLURM scripts (.slurm, .sh):  
Designed for HPC environments using SLURM. Some SLURM jobs call Python or R scripts with the same number in the file name.  
Example:  
05_10_01_BUSCO_all_dinucleotide_calc_parallel_updown_exon_intron_separated.slurm  
→ runs  
05_10_01_BUSCO_all_dinucleotide_calc_parallel_updown_exon_intron_separated.py  
  
Python (.py) and R (.r, .Rmd) scripts:  
Some are standalone, others are executed via SLURM jobs.  
  
⚠️ Note:  
Please adjust directory paths, SLURM resource specifications (cores, memory, time), and environment settings according to your system.  
  
📂 Included Files  
Each analysis directory includes　input files we used for our analyses and several representative output files for reference. 
Some file numbers are missing because certain scripts were removed during code revisions.
  
Author: Yuiki Kondo
