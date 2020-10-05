import haploid.main
import csv
import os
import io
from pathlib import Path

cwd = Path(os.path.realpath(__file__)).parent


def test_main():
    data = Path(cwd / 'data' / 'sample1.csv')
    args = haploid.main.parse_args([str(data)])
    results = haploid.main.run(args)
    haploid.main.write_results(results, args.output_file_path)
    # csvfile = io.StringIO()
    # writer = csv.writer(csvfile)
    # for row in results:
    #     writer.writerow(row)


def test_translate_dna():
    seq = 'CAACAAAAAAAA'
    haploid.main.translateDNA(seq)