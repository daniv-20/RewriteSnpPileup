import argparse
import gzip
import csv

def convert_gzip_to_csv(input_file, output_file):
    """Convert a gzipped file to a CSV file."""
    with gzip.open(input_file, "rt") as gz_file:
        lines = gz_file.readlines()

    with open(output_file, "w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        for line in lines:
            writer.writerow(line.strip().split(","))

    print(f"Data exported to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Convert a gzipped file to a CSV file.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input gzipped file.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output CSV file.")
    args = parser.parse_args()

    convert_gzip_to_csv(args.input, args.output)

if __name__ == "__main__":
    main()
