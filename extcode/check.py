import gzip
import csv

input_file = "/home/nfs/vaithid1/FACETS/RewriteSnpPileup/outputs/TestGzipped.csv.gz"


with gzip.open("/home/nfs/vaithid1/FACETS/RewriteSnpPileup/outputs/TestPythonOut_compArgs.gz", "rt") as f:
    for i, line in enumerate(f):
        print(line.strip())
        if i == 10:  # Stop after the 5th line (index 4 is the 5th line)
            break

