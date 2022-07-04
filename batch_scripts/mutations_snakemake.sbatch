#!/bin/sh
#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20000
#SBATCH --mail-type=END

bash
cd output

#snakemake --forceall --dag -s BEP_pipeline/workflow/Snakefile| dot -Tpdf > dag_proportions.pdf
#snakemake --forceall --rulegraph -s BEP_pipeline/workflow/Snakefile| dot -Tpdf > dag_rulegraph.pdf

snakemake --cores 16 -s ../BEP_pipeline/workflow/Snakefile \
  --resources mem_mb=20000 \
  --configfile /tudelft.net/staff-umbrella/YeastVariantCalling/wsung/BEP_pipeline/config/mutations_config.yaml \
  --config mutations=20 \
  --config experiment_name=mutations_20 --unlock

snakemake --cores 16 -s ../BEP_pipeline/workflow/Snakefile \
  --resources mem_mb=20000 \
  --configfile /tudelft.net/staff-umbrella/YeastVariantCalling/wsung/BEP_pipeline/config/mutations_config.yaml \
  --config mutations=20 \
  --config experiment_name=mutations_20
  
  snakemake --cores 16 -s ../BEP_pipeline/workflow/Snakefile \
  --resources mem_mb=20000 \
  --configfile /tudelft.net/staff-umbrella/YeastVariantCalling/wsung/BEP_pipeline/config/mutations_config.yaml \
  --config mutations=2000 \
  --config experiment_name=mutations_2000 \
  --unlock
 
snakemake --cores 16 -s ../BEP_pipeline/workflow/Snakefile \
  --resources mem_mb=20000 \
  --configfile /tudelft.net/staff-umbrella/YeastVariantCalling/wsung/BEP_pipeline/config/mutations_config.yaml \
  --config mutations=2000 \
  --config experiment_name=mutations_2000 \
  --rerun-incomplete
