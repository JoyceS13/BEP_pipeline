#!/bin/sh
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=30000
#SBATCH --mail-type=END

bash
cd output

#snakemake --forceall --dag -s BEP_pipeline/workflow/Snakefile| dot -Tpdf > dag_proportions.pdf
#snakemake --forceall --rulegraph -s BEP_pipeline/workflow/Snakefile| dot -Tpdf > dag_rulegraph.pdf

snakemake --cores 16 -s ../BEP_pipeline/workflow/Snakefile --resources mem_mb=15000 \
  --configfile /tudelft.net/staff-umbrella/YeastVariantCalling/wsung/BEP_pipeline/config/ploidy5_config.yaml\
  --unlock 

snakemake --cores 16 -s ../BEP_pipeline/workflow/Snakefile --resources mem_mb=30000 \
  --configfile /tudelft.net/staff-umbrella/YeastVariantCalling/wsung/BEP_pipeline/config/ploidy5_config.yaml \
  --rerun-incomplete
