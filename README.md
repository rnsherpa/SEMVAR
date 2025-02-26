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

Input variant list should be in VCF format, ideally with SPDI or rsID in the 3rd column (Check `example/variant_list/example_SNVs.vcf` for a properly formatted example):
- `chrom`: Chromosome (e.g. chr1, chrX, etc.)
- `pos`: 1-based position of variant
- `id`: Unique identifier (usually rsID or SPDI, IGVF Variant-TF Prediction Format requires SPDI)
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
  -v $(pwd)/example/variant_list:/SEMVAR/data/variant_list \
  -v $(pwd)/data/SEMs:/SEMVAR/data/SEMs
  -v path/to/indexed/reference:/SEMVAR/assembly \
  -v $(pwd)/example/example_annotations:/SEMVAR/output \
  semvar -f /SEMVAR/variant_list/example_SNVs.vcf -d /SEMVAR/data/SEMs -a /SEMVAR/assembly/hg38.fa -n 1 -o /SEMVAR/output
```

## Conda example usage

Set up conda environment:

`conda env create -f environment.yml`

Run SEMVAR with example variant list:

`python -m semvar.semvar -f example/variant_list/example_SNVs.vcf -d data/SEMs -a path/to/indexed/reference/hg38.fa -n 1 -o example/example_annotations`