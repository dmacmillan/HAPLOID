#!/home4/dmacmill/python/Python-2.7.2/python
#/Library/Frameworks/Python.framework/Versions/2.7/bin/python

##############################################################
# HAPLOID - Hla Associated PoLymOrphism IDentifier           #
# HAmMOnD - Hla Associated polyMOrphism Discovery            #
#------------------------------------------------------------#
# Author: Daniel MacMillan                                   #
# Contact: drm5@sfu.ca                                       #
##############################################################

import cgi, sys, math

sys.stderr = open("../error-cgi.log", "a")

form = cgi.FieldStorage()
seqs = form.getvalue("sequencesname")
runHAPLOID = form.getvalue("runHAPLOID")

def printHtmlHeaders():
    print "Content-Type: text/html"
    print
    print """<!DOCTYPE html><html><head>
    <script src="//ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script>
    <script src="../cgi-bin/script.js"></script>
    <link rel="stylesheet" href="../css/style.css"></head><body>"""

def printFileHeaders(filename):
    print "Content-Disposition: attachment; filename=\""+filename+"\""
    print "Content-Type:application/octet-stream; name=\""+filename+"\""
    print

codon_dict = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
              'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
              'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
              'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
              'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
              'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
              'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
              'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
              'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
              'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
              'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
              'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
              'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
              'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
              'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
              'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
              '---':'-', 'XXX':'-', '???':'?'}

mixture_dict = {'W':'AT', 'R':'AG', 'K':'GT', 'Y':'CT',
                'S':'CG', 'M':'AC', 'V':'AGC', 'H':'ATC',
                'D':'ATG', 'B':'TGC', 'N':'ATGC', '-':'-'}

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
                nonmix = [x+y for x in nonmix for y in mixture_dict[base]]
        else:
            if (not nonmix):
                nonmix.append(base)
            else:
                nonmix = [x+base for x in nonmix]
    return nonmix

# Flag can be 0, 1, or 2 depending on the desired output
# Flag = 1 will output all mixtures as "X"
# Flag = 2 will output all synonymous mixtures as they are and all non-synonymous mixtures as "X"
# Flag = 3 will output all mixtures in the format [A/B] if a mixture encodes for amino acid A or B
def translateDNA(sequence, resolvecharacter="X", flag=3):
    sequence = sequence.translate(None, ' \n\r\n').upper()
    aaseq = []
    i = 0
    while i < len(sequence):
        try:
            codon = resolveCodon(sequence[i:i+3])
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
                    aaseq.append('['+('/').join(unique)+']')
                else:
                    aaseq.append(unique.pop())
        i += 3
    return aaseq

def parse(inputText):
    val = [x.split('\t') for x in inputText.splitlines()]
    return val

def parseHLA(hla, res=4):
    rval = hla.strip()
    rval = hla.translate(None, ":*")
    rval = hla.upper()
    try:
        int(rval[-1])
    except (ValueError, IndexError) as e:
        rval = rval[:-1]
    return rval[:res+1]

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

def calcMedian(array):
    length = len(array)
    if (length % 2 == 0):
        return ((array[length/2]) + (array[(length/2)-1]))/2
    else:
        return array[length/2]

if (runHAPLOID is not None):
    printHtmlHeaders()
    seqs = getSeqs(parse(seqs))
    print seqs
