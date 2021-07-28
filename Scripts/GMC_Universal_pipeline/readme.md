---
title: "The universal pipelines in Anaconda environment on Gearshift cluster"
author: "Alexander Kurilshikov"
date: "Jul 28, 2021"
---

## Contacts:

This pipeline was created and supported by Ranko Gacesa (rgacesa@gmail.com) and Alex Kurilshikov (alexa.kur@gmail.com). 

## Synopsis:

The idea of this pipeline is to create universal, controlled system within GMC to run our standard things for metagenomic data processing (QC, taxonomy and functional binning, other pipelines) is a simple and straightforward way on **Gearshift** cluster. The system also can be transferred to other Servers if necessary. 

## What is implemented now:

1. Biobakery tools (Kneaddata, Metaphlan 3, HUMANN 3, ShortBRED)

## What you need:

1. Access to **Gearshift** cluster
2. Access to **umcg-tifn** group on the cluster

## Stuff location:

All code and databases are located here:

```
/groups/umcg-tifn/tmp01/tools
```

## Before you start

There are only one thing you need to do **once** to connect to our collaborative *Conda* environment. Just run:

```
$ bash /groups/umcg-tifn/tmp01/tools/activate_biobakery_conda.sh
```

after that, you just need to check that:
1. your personal  ```~/.condarc``` file, looks like this:

```
$ cat ~/.condarc

channels:
  - biobakery
  - conda-forge
  - bioconda
  - defaults
envs_dirs:
  - /groups/umcg-tifn/tmp01/tools/conda_microbiome_July2021/envs
pkgs_dirs:
  - /groups/umcg-tifn/tmp01/tools/conda_microbiome_July2021/pkgs
```

2. your personal ```~/.bashrc``` file has one specific line:

```
$ cat ~/.bashrc |grep conda

. /apps/software/Anaconda3/5.3.0/etc/profile.d/conda.sh
```

If everything goes fine, you can logout and login again, and your *Conda* environment will be set up properly. When you use something from *Conda*, just don't forget to load corresponding module and activate the proper environment. Also you should do the same if you run the script via **SLURM** (i.e. using *sbatch* command)

```
$ module load Anaconda3/5.3.0
$ conda activate tifn_biobakery_July2021
```

Currently, only **Biobakery** environment is available, but that will change in the future.

## Using biobakery tools in the environment

Here comes Ranko part
