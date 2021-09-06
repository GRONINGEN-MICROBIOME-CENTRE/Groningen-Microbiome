Title: Structural variation calling for metagenomic data

Author: Daoming Wang

E-mail: wangdaoming94@outlook.com

Date: Nov 2019

# 1. Abstract

Just like human genome, there are also lots of different types of genetic variants in microbial genomes, including SNPs, CNVs, and SVs etc, which contribute to the phenotypic diversity of microbial strains. Many bioinformatic tools have been devised to digging these variants, of which, SGVFinder is a tool which enabled the discovery of SVs from shot-gun metagenomic sequencing data of environmental microbiome. Two classes have been defined for the metagenomic SVs: deletion SVs and variable SVs. SGVFinder will not work on a single sample, and large cohort-level samples is better for it.

# 2. Preparation

## 2.1 Requirements

- **GEM Mapper**, this version will not work with GEM3, so please try the older versions, for instance version 2:

```shell
wget https://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%202/GEM-binaries-Linux-x86_64-core_2-20121106-022124.tbz2
bzip2 -dc GEM-binaries-Linux-x86_64-core_2-20121106-022124.tbz2 | tar -xvf -
```

The binaries should be added into your PATH variables, so that a simple call to "gem-mapper" would be successful, add the codes below into ~/.bashrc:

```shell
export PATH="your_path/GEM-binaries-Linux-x86_64-core_2-20121106-022124:$PATH"
 
Then source your ~/.bashrc file to update the environment variables:
source ~/.bashrc
```

- **Python     2.7.8**     and following modules:

- - numpy (tested with 1.14.2)
  - biopython (tested with 1.68)
  - ujson (tested with 1.35)
  - pandas (tested with 0.23.4)
  - scipy (tested with 1.1.0)
  - bokeh (tested with 0.12.6)

Strongly recommend using Anaconda, it can provide a user-friendly environment for python module management and language version switch.

- **Cpp11**
- **Cython**, it can also been installed with Anaconda

## 2.2 Download and Install

Download the tools and databases, and install the tools:

```shell
wget https://zenodo.org/record/3237975/files/DataFiles.tar.gz
gzip -dc DataFiles.tar.gz > DataFiles
git clone https://github.com/segalab/SGVFinder.git
mv DataFiles SGVFinder/
cd SGVFinder/src/cy_ext/
python setup.py build_ext
```

**Note:**

There are several bugs need to be corrected, otherwise errors will be generated when running this tool.

```python
(1) SGVFinder.py 
line 32: ind1 = int((pos1 + (average_read_length / 2)) / bin_size) -> ind1 = int((int(pos1) + (int(average_read_length) / 2)) / bin_size)
line 35: ind2 = int((pos2 + (average_read_length / 2)) / bin_size) -> ind2 = int((int(pos2) + (int(average_read_length) / 2)) / bin_size)
 
(2) SGVF_cmd.py
line 39: if args.byother: -> if args.byorig:
line 15: 
parser.add_argument('--x_coverage', help = 'The desired coverage across the genome in units of 100bp reads. This parameter is used to determine bin size: bin_size = rate_param/x_coverage (Default = 0.1)', type=float, default = 0.1) -> 
parser.add_argument('--x_coverage', help = 'The desired coverage across the genome in units of 100bp reads. This parameter is used to determine bin size: bin_size = rate_param/x_coverage (Default = 0.01)', type=float, default = 0.01) # this is not a bug, but this default parameter isn't equal with the corresponding default parameter setting in SGVF_PerFile.py, which may misleading user, then cause error (users tend to use default parameters), you can also set --x_coverage=0.01 to avoid error when running SGVF_cmd.py.
```

 

# 3. Input Files

(1) Single end reads:

- input_a.fastq
- input_b.fastq
- input_c.fastq

 

(2) Paired end reads:

- input_a_1.fastq,     input_a_2.fastq
- input_b_1.fastq,     input_b_2.fastq
- input_c_1.fastq,     input_c_2.fastq

 

# 4. Call SVs

## 4.1 Run ICRA

Fine map the metagenomic reads to avoid ambiguous alignment.

For single end reads:

```shell
mkdir -p tmp/icra
python ICRA_cmd.py tmp/icra input_a 
```

For paired end reads:

```shell
mkdir -p tmp/icra
python ICRA_cmd.py tmp/icra input_a --pe
```

*For this step, peak memory 22 Gb per file, run time 7 hrs per file.*

**Note:**

- Only accept unzipped fastQ file;
- 'tmp/icra' is a directory, not a file.

## 4.2 Run SGVF_PerFile

Generates coverage maps per bacteria per sample, it will generate a automatically from a glob string with the command-line wrapper:

```shell
mkdir -p tmp/SGVF_PerFile
python SGVF_PerFile_cmd.py tmp/icra/input_a.jsdel tmp/SGVF_PerFile/input.map.jsdel 100 --x_coverage 0.01 --rate_param 10
```

*For this step, peak memory 10 Gb per run (1 sample), run time 0.1 hrs per file*.

## 4.3 Run SGVF

**After all input fastq files including input_a, input_b, input_c……have been processed** by above steps, all of their coverage map results should in the directory 'tmp/coverage_map',  SGVFinder **WILL NOT** work on a single sample.

```shell
mkdir -p results
mkdir -p results/html
python SGVF_cmd.py input_glob_string results/dsgv.csv results/vsgv.csv --min_samp_cutoff 5 --x_coverage 0.01 --rate_param 10 --browser_path results/html --csv_output
```

*For this step, peak memory 5 Gb per run (10 samples), run time 0.5 hrs per run (10 samples).*

**Note:**

- For small cohort, reference-based method is recommended, add parameter --byorig ;

- input_glob_string is a quoted File name regular expression, such as: "/my/path/*.jsdel";

- --min_samp_cutoff also should be changed to suitable with your population, 10% * population sample size is recommanded.

  

# 5. Final results

- results/dsgv.csv: deletion SVs profile
- results/vsgv.csv: variable SVs profile

 

# 6. After that…

## 6.1 Transformation

If you didn't use --csv_output, you can transform the file format from pandas dataframe to CSV with the script below:

```python
# python
import pandas as pd
dsgv_df = pd.read_pickle("results/dsgv.df")
dsgv_df.to_csv("results/dsgv.csv")
vsgv_df = pd.read_pickle("results/vsgv.df")
vsgv_df.to_csv("results/vsgv.csv")
```

## 6.2 Results description

- results/dsgv/dsgv.csv: Binomial type, 1/0 values represent the status of present/absent for deletion SVs
- results/vsgv/vsgv.csv: Ratio type, continues values represent the coverage of variable SVs

 

# 7. Reference

- Zeevi D, Korem T, Godneva A, Bar N, Kurilshikov A, Lotan-Pompan M, Weinberger A, Fu J, Wijmenga C, Zhernakova A & Segal E (2019) Structural variation in the gut microbiome associates with host health. Nature 568, 43–48.
- ICRA and SGVFinder: https://github.com/segalab/SGVFinder/
- Exmple in python script: https://github.com/segalab/SGVFinder/blob/master/src/linear_example.py
