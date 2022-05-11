ml Anaconda3/5.3.0
#source activate /groups/umcg-lld/tmp01/other-users/umcg-sandreusanchez/Binning_HIV/Binning_env
source activate /groups/umcg-lld/tmp01/other-users/umcg-sandreusanchez/Binning_HIV/vamp
snakemake --jobs 20 --configfile config.json --snakefile /groups/umcg-lld/tmp01/other-users/umcg-sandreusanchez/Binning_HIV/workflow/vamb.snake.conda.py --latency-wait 60 --use-conda --cluster 'sbatch -t {cluster.walltime} --mem={cluster.mem} --cpus-per-task={cluster.ppn} --error={cluster.error} --job-name={cluster.name} --output={cluster.output}'  --cluster-config cluster.json  --cluster-status 'python slurm-status.py' --conda-base-path /groups/umcg-lld/tmp01/other-users/umcg-sandreusanchez/Binning_HIV/Binning_env/bin/ #--conda-frontend mamba

# --cluster 'sbatch -t {cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.c} --error={cluster.error} --job-name={cluster.name} --output={cluster.output}'  --cluster-config cluster.json --snakefile Analyze_samples.smk --cluster-status 'python slurm-status.py'


#"qsub -l time={params.walltime} -l cpu-per-task={params.ppn} -l mem={params.mem}" #'sbatch -t {params.walltime} --mem={params.mem} --nodes={params.ppn}'

#"qsub -l walltime={params.walltime} -l nodes=1:ppn={params.ppn} -l mem={params.mem}"
#'sbatch -t {params.walltime} --mem={params.mem} --nodes=1:ppn={params.ppn}'

