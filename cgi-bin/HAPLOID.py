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
patients = form.getvalue("sequencesname")
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
    rval = hla.strip().translate(None, "*:").upper()
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

def buildUniqueHlas(patients):
    hlas = set()
    for patient in patients:
        [hlas.add(x) for x in patients[patient]['A']]
        [hlas.add(x) for x in patients[patient]['B']]
        [hlas.add(x) for x in patients[patient]['C']]
    return list(hlas)
    
def buildHLAHash(uniquehlas):
    res = {}
    for i,hla in enumerate(uniquehlas):
        res[hla] = i
    return res

def analyzeHLAs(patients, uniquehlas):
    results = {}
    for uhla in uniquehlas:
        results[uhla] = {}
        for patient in patients:
            n = len(patients)
            for pos, aa in enumerate(patients[patient]['seq']):
                #if (len(aa) > 1):
                    #continue
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
    return results
    
def displayResults(analysis, num):
    #print '{}<br>'.format(analysis)
    print '''<table id="output_table">
            <th>HLA</th>
            <th>CODON</th>
            <th>AA</th>
            <th>DIRECTION</th>
            <th>TT</th>
            <th>TF</th>
            <th>FT</th>
            <th>FF</th>
            <th>N</th>
            <th>ODDS RATIO</th>
            <th>P-VALUE</th>'''
    for uhla in analysis:
        for pos in analysis[uhla]:
            num_patients = num
            #print 'looking at position: {}<br>'.format(pos)
            aas = set([x for x in analysis[uhla][pos]['tt']]+[x for x in analysis[uhla][pos]['ft']])
            if (len(aas) == 1):
                num_patients -= 1
                continue
            for aa in aas:
                if (aa[0] == '['):
                    continue
                #print '&nbsp;looking at: {}<br>'.format(aa)
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
                tf = sum([analysis[uhla][pos]['tt'][x] for x in analysis[uhla][pos]['tt'] if x != aa])
                if aa in mixtt:
                    tf += mixtt[aa]
                if (aa in analysis[uhla][pos]['ft']):
                    ft = analysis[uhla][pos]['ft'][aa]
                else:
                    ft = 0
                ff = sum([analysis[uhla][pos]['ft'][x] for x in analysis[uhla][pos]['ft'] if x != aa])
                if aa in mixft:
                    ff -= mixft[aa]
                    num_patients -= 1
                try:
                    OR = (float(tt*ff)/(tf*ft))
                except ZeroDivisionError:
                    OR = 'Indeterminate'
                try:
                    abfact = math.factorial(tt+tf)
                    cdfact = math.factorial(ft+ff)
                    acfact = math.factorial(tt+ft)
                    bdfact = math.factorial(tf+ff)
                    afact = math.factorial(tt)
                    bfact = math.factorial(tf)
                    cfact = math.factorial(ft)
                    dfact = math.factorial(ff)
                    nfact = math.factorial(tt+tf+ft+ff)
                    logp = math.log(abfact) + math.log(cdfact) + math.log(acfact) + math.log(bdfact) - math.log((afact * bfact * cfact * dfact * nfact))
                    pval = round(math.exp(logp),8)
                except ZeroDivisionError:
                    pval = "Indeterminate"
                if (OR < 1):
                    direction = 'non-adapted'
                elif (OR == 'indeterminate') or (OR > 1):
                    direction = 'adapted'
                else:
                    direction = 'n/a'
                print '<tr>'
                print '<td>{}</td>'.format(uhla)
                print '<td>{}</td>'.format(pos+1)
                print '<td>{}</td>'.format(aa)
                print '<td>{}</td>'.format(direction)
                print '<td>{}</td>'.format(tt)
                print '<td>{}</td>'.format(tf)
                print '<td>{}</td>'.format(ft)
                print '<td>{}</td>'.format(ff)
                print '<td>{}</td>'.format(num_patients)
                print '<td>{}</td>'.format(OR)
                print '<td>{}</td>'.format(pval)
                print '</tr>'
    print '</table>'

if (runHAPLOID is not None):
    printHtmlHeaders()
    patients = getSeqs(parse(patients))
    #print patients
    uniquehlas = buildUniqueHlas(patients)
    #print uniquehlas
    results = analyzeHLAs(patients, uniquehlas)
    displayResults(results,len(patients))