#!/bin/bash
#
#$ -S /bin/bash  # defines bash as the shell for execution
#$ -N cross      # name of the command that will be listed in the queue
#$ -cwd          # change to the current directory
#$ -j y          # join error and standard output in one file
#$ -o ~/log      # standard output file or directory
#$ -q archgen.q  # queue
#$ -pe smp 2     # use X CPU cores
#$ -l h_vmem=5G  # request XGb of memory
#$ -V            # load personal profile
#$ -t 1-225      # array job length
#$ -tc 25        # number of concurrently submitted tasks

date

echo Task in Array: ${SGE_TASK_ID}
i=$((SGE_TASK_ID - 1))

ds_to_explore=($(seq 100 100 1500))
dt_to_explore=($(seq 100 100 1500))
dss=()
dts=()
for ds in "${ds_to_explore[@]}"
do
  for dt in "${dt_to_explore[@]}"
  do
    dss+=($ds)
    dts+=($dt)
  done
done

current_ds=${dss[${i}]}
current_dt=${dts[${i}]}

echo ds: ${current_ds}
echo dt: ${current_dt}

apptainer exec \
  --bind=/path/to/your/analysis/directory \
  path/to/your/apptainer_mobest.sif \
  Rscript path/to/your/mobestRscript.R ${i} ${current_ds} ${current_dt} \
  /

date
exit 0
