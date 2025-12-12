cellranger-arc + cellbender pipeline.
requires folders for each sample (rna+atac library pair) with a `libraries.csv` file present in each.

the launcher script uses a `samples.txt` file found in the same folder. it
just contains one line with the folder name (each of which contains the
`libraries.csv` file).
