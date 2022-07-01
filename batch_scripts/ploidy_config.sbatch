#!/bin/sh
#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=15000
#SBATCH --mail-type=END

bash
cd output

#snakemake --forceall --dag -s BEP_pipeline/workflow/Snakefile| dot -Tpdf > dag_proportions.pdf
#snakemake --forceall --rulegraph -s BEP_pipeline/workflow/Snakefile| dot -Tpdf > dag_rulegraph.pdf

snakemake --cores 16 -s BEP_pipeline/workflow/Snakefile --resources mem_mb=15000 \
  --configfile /tudelft.net/staff-umbrella/YeastVariantCalling/wsung/BEP_pipeline/config/ploidy3_config.yaml
  
#snakemake --cores 16 -s BEP_pipeline/workflow/Snakefile --resources mem_mb=15000 \
#  --configfile /tudelft.net/staff-umbrella/YeastVariantCalling/wsung/BEP_pipeline/config/ploidy4_config.yaml

#snakemake --cores 16 -s BEP_pipeline/workflow/Snakefile --resources mem_mb=15000 \
#  --configfile /tudelft.net/staff-umbrella/YeastVariantCalling/wsung/BEP_pipeline/config/ploidy5_config.yaml
