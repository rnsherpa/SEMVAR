# SEMVAR: SEM Variant Annotator

## Usage

'''
usage: SEM_variant_annotation [-h] [--file FILE] [--semsdir SEMSDIR]
                              [--assembly ASSEMBLY] [--baselines BASELINES]
                              [--n_processes N_PROCESSES] [--outdir OUTDIR]
'''

| Argument | Description |
| -------- | ----------- |
| --file   | Input variant list |
| --semsdir | Path to directory containing SEMs |
| --assembly | Path to indexed reference genome |
| --baselines | Path to file containing baseline values for each SEM |
| --n_processes | Number of processes to use |
| --outdir | Output path for variant annotations |

## Example

Set up conda environment:

`conda env create -f environment.yml`

Run SEMVAR with example variant list:

`python SEMVAR.py -f example/variant_list/example_SNVs.tsv -s data/SEMs/ -a path/to/indexed/reference/hg38.fa -b data/BASELINE/SEMs_baseline_norm.txt -n 1 -o example/example_annotations`

