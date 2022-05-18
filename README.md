# nf-compare-miRNome

Nextflow pipeline that compare miRNA targets of miRNA and 3'UTR reference and mutate sequences.

Basic idea:

1)  Create a consensus FASTA sequence from VCF file and extract the reference and mutated miRNA and 3'UTR sequences

2)  Use miRmap and TargetScan to predict miRNA targets of reference and mutate sequences

3)Compare microRNA targets and calcule changes to plot microRNA target changes

------------------------------------------------------------------------

### Workflow overview

![General Workflow](dev_notes/Workflow.png)

------------------------------------------------------------------------

## Requirements

#### Compatible OS\*:

-   [Ubuntu 20.04]

#### Software:

|                                     Requirement                                     |          Version           |           Required Commands \*            |
|:-----------------------------------------------------------------------------------:|:--------------------------:|:-----------------------------------------:|
|                  [Plan9 port](https://github.com/9fans/plan9port)                   | Latest (as of 10/01/2019 ) |                  mk \*\*                  |
|                        [Nextflow](https://www.nextflow.io/)                         |          21.04.2           |                   run\*                   |
|                           [R](https://www.r-project.org/)                           |           3.4.4            |               See R scripts               |
|                 [Python](https://www.python.org/downloads/source/)                  |           3.6.9            |            See python scripts             |
|                                      bcftools                                       |           1.10.2           |               view, index\*               |
|                                      bedtools                                       |           2.27.1           |                getfasta\*                 |
|                                        bgzip                                        |          1.10.2-3          |                  bgzip\*                  |
|                 [Python](https://www.python.org/downloads/source/)                  |           3.6.9            |        \*\*\*See mirmap_script.py         |
|                     [miRmap](https://pypi.org/project/mirmap/)                      |           0.0.1            |        \*\*\*See mirmap_script.py         |
|                  [Biopython](https://biopython.org/wiki/Download)                   |            1.79            |    \*\*\*SeqIO (See mirmap_script.py)     |
| [TargetScan](http://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi) |            7.0             | \*Script included in the pipeline modules |

\* These commands must be accessible from your `$PATH` (*i.e.* you should be able to invoke them from your command line).

\*\* Plan9 port builds many binaries, but you ONLY need the `mk` utility to be accessible from your command line.

### Installation

Download compare-miRNA-pairs.nf from Github repository:

    git clone https://github.com/Ed-G655/nf-compare-miRNome.git

------------------------------------------------------------------------

#### Test

To test compare-miRNA-pairs.nf execution using test data, run:

    bash runtest.sh

Your console should print the Nextflow log for the run, once every process has been submitted, the following message will appear:

     ======
     Basic pipeline TEST SUCCESSFUL
     ======

compare-miRNA-pairs.nf results for test data should be in the following file:

    nf-compare-miRNome/test/results/

------------------------------------------------------------------------

### Usage

To run compare-miRNA-pairs go to the pipeline directory and execute:

    nextflow run nf-compare-miRNome.nf --mirnabed <path to input 1> --utrbed <path to input 2>
      --vcf <path to input 3> --fasta <path to input >   [--output_dir path to results ]

For information about options and parameters, run:

    nextflow run nf-compare-miRNome.nf --help

------------------------------------------------------------------------

#### References

#### Autors

José Eduardo García López
