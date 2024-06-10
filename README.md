# SEMVAR: SEM Variant Annotator

## Usage

```
usage: python semvar.py [-h] [--file FILE] [--semsdir SEMSDIR]
                              [--assembly ASSEMBLY] [--baselines BASELINES]
                              [--n_processes N_PROCESSES] [--outdir OUTDIR]
```

| Argument | Description |
| -------- | ----------- |
| --file   | Input variant list |
| --semsdir | Path to directory containing SEMs |
| --assembly | Path to indexed reference genome |
| --baselines | Path to file containing baseline values for each SEM |
| --n_processes | Number of processes to use for Python multiprocessing |
| --outdir | Output path for variant annotations |

Input variant list should have 5 columns (Check `example/variant_list/example_SNVs.tsv` for a properly formatted example):
- `chrom`: Chromosome (e.g. chr1, chrX, etc.)
- `start`: 0-based start position of variant
- `end`: 0-based end position of variant
- `ref`: Reference allele
- `alt`: Alternate allele

## Docker example usage (Recommended)

Build docker image (Temporary until image hosted):

`docker build -t semvar .`

Run SEMVAR docker image:
```
docker run -it --rm \
  -e HOST_UID=$(id -u) \
  -e HOST_GID=$(id -g) \
  -v $(pwd)/example/variant_list:/SEMVAR/variant_list \
  -v path/to/indexed/reference:/SEMVAR/assembly \
  -v $(pwd)/example/example_annotations:/SEMVAR/output \
  semvar -f /SEMVAR/variant_list/example_SNVs.tsv -s /SEMVAR/data/SEMs -a /SEMVAR/assembly/hg38.fa -b /SEMVAR/data/BASELINE/SEMs_baseline_norm.txt -n 1 -o /SEMVAR/output
```

## Conda example usage

Set up conda environment:

`conda env create -f environment.yml`

Run SEMVAR with example variant list:

`python SEMVAR.py -f example/variant_list/example_SNVs.tsv -s data/SEMs -a path/to/indexed/reference/hg38.fa -b data/BASELINE/SEMs_baseline_norm.txt -n 1 -o example/example_annotations`