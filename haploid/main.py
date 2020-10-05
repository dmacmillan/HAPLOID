import argparse
import sys
import csv
import os

from pathlib import Path
from scipy import stats

if __name__ == '__main__':
    main()

codon_dict = {
    'TTT': 'F',
    'TTC': 'F',
    'TTA': 'L',
    'TTG': 'L',
    'TCT': 'S',
    'TCC': 'S',
    'TCA': 'S',
    'TCG': 'S',
    'TAT': 'Y',
    'TAC': 'Y',
    'TAA': '*',
    'TAG': '*',
    'TGT': 'C',
    'TGC': 'C',
    'TGA': '*',
    'TGG': 'W',
    'CTT': 'L',
    'CTC': 'L',
    'CTA': 'L',
    'CTG': 'L',
    'CCT': 'P',
    'CCC': 'P',
    'CCA': 'P',
    'CCG': 'P',
    'CAT': 'H',
    'CAC': 'H',
    'CAA': 'Q',
    'CAG': 'Q',
    'CGT': 'R',
    'CGC': 'R',
    'CGA': 'R',
    'CGG': 'R',
    'ATT': 'I',
    'ATC': 'I',
    'ATA': 'I',
    'ATG': 'M',
    'ACT': 'T',
    'ACC': 'T',
    'ACA': 'T',
    'ACG': 'T',
    'AAT': 'N',
    'AAC': 'N',
    'AAA': 'K',
    'AAG': 'K',
    'AGT': 'S',
    'AGC': 'S',
    'AGA': 'R',
    'AGG': 'R',
    'GTT': 'V',
    'GTC': 'V',
    'GTA': 'V',
    'GTG': 'V',
    'GCT': 'A',
    'GCC': 'A',
    'GCA': 'A',
    'GCG': 'A',
    'GAT': 'D',
    'GAC': 'D',
    'GAA': 'E',
    'GAG': 'E',
    'GGT': 'G',
    'GGC': 'G',
    'GGA': 'G',
    'GGG': 'G',
    '---': '-',
    'XXX': '-',
    '???': '?'
}

mixture_dict = {
    'W': 'AT',
    'R': 'AG',
    'K': 'GT',
    'Y': 'CT',
    'S': 'CG',
    'M': 'AC',
    'V': 'AGC',
    'H': 'ATC',
    'D': 'ATG',
    'B': 'TGC',
    'N': 'ATGC',
    '-': '-'
}


def resolveCodon(codon):
    nonmix = []
    if (codon in codon_dict):
        return [codon]
    elif (codon.count('-') + codon.count('X') == 3):
        return ['---']
    elif (1 <= codon.count('-') <= 2) or (1 <= codon.count('X') <= 2):
        return ['???']
    for base in codon:
        # Check for mixtures
        if (base in mixture_dict):
            if (not nonmix):
                nonmix = [x for x in mixture_dict[base]]
            else:
                nonmix = [x + y for x in nonmix for y in mixture_dict[base]]
        else:
            if (not nonmix):
                nonmix.append(base)
            else:
                nonmix = [x + base for x in nonmix]
    return nonmix


# Flag can be 0, 1, or 2 depending on the desired output
# Flag = 1 will output all mixtures as "X"
# Flag = 2 will output all synonymous mixtures as they are and all non-synonymous mixtures as "X"
# Flag = 3 will output all mixtures in the format [A/B] if a mixture encodes for amino acid A or B
def translateDNA(sequence, resolvecharacter="X", flag=3):
    chars = [' ', '\r\n', '\n', '\r']
    for char in chars:
        if char in sequence:
            sequence = sequence.replace(char, '')
    sequence = sequence.upper()
    aaseq = []
    i = 0
    while i < len(sequence):
        try:
            codon = resolveCodon(sequence[i:i + 3])
        except IndexError:
            codon = resolveCodon('???')
        # If the codon has no mixture bases just add it to the amino acid chain
        if len(codon) <= 1:
            aaseq.append(codon_dict[codon[0]])
        # Codon contains mixture base
        else:
            # If flag is set to 1
            if (flag == 1):
                aaseq.append(resolvecharacter)
            # If flag is set to 2
            elif (flag == 2):
                unique = set([codon_dict[potential] for potential in codon])
                # If there is more than resolved one amino acid
                if (len(unique) > 1):
                    aaseq.append(resolvecharacter)
                else:
                    aaseq.append(unique.pop())
            # If flag is set to 3
            else:
                unique = set([codon_dict[potential] for potential in codon])
                # If there is more than resolved one amino acid
                if (len(unique) > 1):
                    aaseq.append('[' + ('/').join(unique) + ']')
                else:
                    aaseq.append(unique.pop())
        i += 3
    return aaseq


def parse(inputText):
    val = [x.split('\t') for x in inputText.splitlines()]
    return val


def parseCSV(inputText):
    val = [x.split(',') for x in inputText.splitlines()]
    return val


def parseHLA(hla, res=4):
    chars = ['*', ':']
    for char in chars:
        if char in hla:
            hla = hla.replace(char, '')
    rval = hla.strip().upper()
    try:
        int(rval[-1])
    except (ValueError, IndexError) as e:
        rval = rval[:-1]
    return rval[:res + 1]


def getSeqs(seqs):
    d = {}
    for patient in seqs:
        pid = patient[0]
        nuseq = patient[-1]
        aaseq = translateDNA(nuseq)
        d[pid] = {'A': [], 'B': [], 'C': [], 'seq': aaseq}
        for hla in patient[1:-1]:
            hla = parseHLA(hla)
            if (hla == ""):
                continue
            if (hla[0].upper() == 'A'):
                d[pid]['A'].append(hla)
            elif (hla[0].upper() == 'B'):
                d[pid]['B'].append(hla)
            else:
                d[pid]['C'].append(hla)
    return d


def buildUniqueHlas(patients):
    hlas = set()
    for patient in patients:
        [hlas.add(x) for x in patients[patient]['A']]
        [hlas.add(x) for x in patients[patient]['B']]
        [hlas.add(x) for x in patients[patient]['C']]
    return list(hlas)


def analyzeHLAs(patients, uniquehlas, hlacountfilter=None, aacountfilter=1):
    results = {}
    seqcounts = {}
    for patient in patients:
        for pos, aa in enumerate(patients[patient]['seq']):
            if pos not in seqcounts:
                seqcounts[pos] = {}
            if aa in seqcounts[pos]:
                seqcounts[pos][aa] += 1
            else:
                seqcounts[pos][aa] = 1
    for uhla in uniquehlas:
        hlacount = 0
        results[uhla] = {}
        for patient in patients:
            n = len(patients)
            if uhla in patients[patient][uhla[0]]:
                hlacount += 1
            for pos, aa in enumerate(patients[patient]['seq']):
                if (seqcounts[pos][aa] < aacountfilter):
                    continue
                #if (len(aa) > 1):
                #continue
                if (aa in ['X', '-']):
                    continue
                if pos not in results[uhla]:
                    results[uhla][pos] = {'tt': {}, 'ft': {}, 'n': n}
                if uhla not in patients[patient][uhla[0]]:
                    if aa not in results[uhla][pos]['ft']:
                        results[uhla][pos]['ft'][aa] = 1
                    else:
                        results[uhla][pos]['ft'][aa] += 1
                else:
                    if aa not in results[uhla][pos]['tt']:
                        results[uhla][pos]['tt'][aa] = 1
                    else:
                        results[uhla][pos]['tt'][aa] += 1
        if ((hlacountfilter) and
            (hlacount < hlacountfilter)) or (hlacount == len(patients)):
            #print 'There are only {} patient/s with hla {} which is < {}, removing.<br>'.format(hlacount,uhla,countfilter)
            results.pop(uhla, None)
    return results


def getResults(analysis):
    results = []
    for uhla in analysis:
        for pos in analysis[uhla]:
            aas = set([x for x in analysis[uhla][pos]['tt']] +
                      [x for x in analysis[uhla][pos]['ft']])
            copy = list(aas)
            for aa in copy:
                if (aa[0] == '['):
                    copy.remove(aa)
            if (len(aas) == 1) or (len(copy) == 1):
                continue
            for aa in aas:
                if (aa[0] == '['):
                    continue
                mixft = {}
                for pot in analysis[uhla][pos]['ft']:
                    if pot[0] == '[':
                        t = pot[1:-1].split('/')
                        for i in t:
                            mixft[i] = analysis[uhla][pos]['ft'][pot]
                mixtt = {}
                for pot in analysis[uhla][pos]['tt']:
                    if pot[0] == '[':
                        t = pot[1:-1].split('/')
                        for i in t:
                            mixtt[i] = analysis[uhla][pos]['tt'][pot]
                if (aa in analysis[uhla][pos]['tt']):
                    tt = analysis[uhla][pos]['tt'][aa]
                else:
                    tt = 0
                tf = sum([
                    analysis[uhla][pos]['tt'][x]
                    for x in analysis[uhla][pos]['tt'] if x != aa
                ])
                if aa in mixtt:
                    tf -= mixtt[aa]
                if (aa in analysis[uhla][pos]['ft']):
                    ft = analysis[uhla][pos]['ft'][aa]
                else:
                    ft = 0
                ff = sum([
                    analysis[uhla][pos]['ft'][x]
                    for x in analysis[uhla][pos]['ft'] if x != aa
                ])
                if aa in mixft:
                    ff -= mixft[aa]
                n = tt + tf + ft + ff
                OR, pval = stats.fisher_exact([[tt, tf], [ft, ff]])
                if (OR < 1):
                    direction = 'non-adapted'
                elif (OR == 'indeterminate') or (OR > 1):
                    direction = 'adapted'
                else:
                    direction = 'n/a'
                results.append([
                    uhla,
                    str(pos + 1), aa, direction,
                    str(tt),
                    str(tf),
                    str(ft),
                    str(ff),
                    str(n),
                    str(OR),
                    str(pval)
                ])
    return results


def run(args):
    text = None
    with open(args.data_csv) as f:
        text = f.read()
    patients = getSeqs(parseCSV(text))
    unique_hlas = buildUniqueHlas(patients)
    results = analyzeHLAs(patients, unique_hlas, args.hla_count_filter,
                          args.amino_acid_count_filter)
    results = getResults(results)
    results = sorted(results, key=lambda x: (x[2]))
    results = sorted(results, key=lambda x: (x[3]), reverse=True)
    results = sorted(results, key=lambda x: (int(x[1]), x[0]))
    return results


def write_results(results, outfile):
    fieldnames = [
        'HLA',
        'Codon',
        'Amino Acid',
        'Direction',
        'TT',
        'TF',
        'FT',
        'FF',
        'n',
        'Odds Ratio',
        'Pvalue',
        'Qvalue',
    ]
    print(f'Writing results to: "{outfile}"')
    with open(outfile, 'w', newline='') as o:
        writer = csv.writer(o)
        writer.writerow(fieldnames)
        for row in results:
            writer.writerow(row)


def main():
    args = parse_args(sys.argv[1:])
    results = run(args)
    write_results(results, args.output_file_path)


def parse_args(args):
    parser = argparse.ArgumentParser(
        description='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'data_csv',
        help=
        'Do not include headers. Format is one patient per line, and exactly 8 columns per patient. Each patient should have their unique patient ID in the first column. The following 6 columns should be HLA types. The final column must be the nucleotide sequence of the patient\'s HIV gene of interest.'
    )
    parser.add_argument(
        '-o',
        '--output_file_path',
        type=Path,
        default=(Path(os.getcwd()) / 'results.csv').resolve(),
        help=
        'The path to the file to write results to. Directory structure must already exist!'
    )
    parser.add_argument(
        '--hla_count_filter',
        type=int,
        default=1,
        help=
        'Exclude HLAs from analysis if number of patients with said HLAs is less than this number. Minimum = 1, required.'
    )
    parser.add_argument(
        '--amino_acid_count_filter',
        type=int,
        default=1,
        help=
        'Exclude amino acids (AA)s from analysis if number of patients with said AA at any given position is less than this number. Minimum = 1, required.'
    )
    parser.add_argument(
        '--max_pvalue',
        type=float,
        default=0.05,
        help='Do not output p-values larger than this value. Default = 0.05.')
    parser.add_argument(
        '--max_qvalue',
        type=float,
        default=0.02,
        help='Do not output q-values larger than this value. Default = 0.05.')
    return parser.parse_args(args)
