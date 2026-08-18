"""Q0-B evidence: BRAF historical numbering 3-nt difference.
Fetches NP_004324.1 (historical RefSeq) and M95712.1 CDS translation,
aligns against the current canonical P15056, records the差异.
"""
import json, re, urllib.request
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent

UA = {'User-Agent': 'MetaSieve-StageX0c/1.0 (research audit)'}

def http_get(url, tries=3):
    import time
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode('utf-8')
        except Exception as e:
            last = e
            time.sleep(2)
    raise last

def parse_fasta(text):
    seq = ''.join(l.strip() for l in text.splitlines() if not l.startswith('>'))
    header = next((l for l in text.splitlines() if l.startswith('>')), '')
    return header, seq

def translate(seq):
    code = {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
    out = []
    for i in range(0, len(seq) - 2, 3):
        out.append(code.get(seq[i:i+3], 'X'))
    return ''.join(out)

# current canonical
p15056 = ''.join(l.strip() for l in (PARENT/'uniprot'/'P15056.fasta').read_text().splitlines() if not l.startswith('>'))

# historical RefSeq NP_004324.1
h1, s1 = parse_fasta(http_get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=NP_004324.1&rettype=fasta&retmode=text'))
# 1992 mRNA M95712.1 CDS
h2, cdna = parse_fasta(http_get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=M95712.1&rettype=fasta_cds_na&retmode=text'))
t2 = translate(cdna).rstrip('*')

print('header NP_004324.1:', h1[:90])
print('len NP_004324.1:', len(s1), 'len M95712 CDS translation:', len(t2), 'len P15056:', len(p15056))

def diff_positions(a, b):
    n = min(len(a), len(b))
    return [i for i in range(n) if a[i] != b[i]]

d1 = diff_positions(s1, p15056)
d2 = diff_positions(t2, p15056)
print('NP_004324.1 vs P15056 differing positions (aa):', d1[:20], 'count', len(d1))
print('M95712 translation vs P15056 differing positions:', d2[:20], 'count', len(d2))
for i in (d1 or [])[:6]:
    lo, hi = max(0, i-6), min(len(p15056), i+8)
    print(f'  pos {i+1}: hist={s1[lo:hi]} canon={p15056[lo:hi]}')

out = {
    'schema': 'MetaSieve.StageX0c.BRAFHistoricalEvidence.v1',
    'braf_current_canonical': 'UniProt P15056 SV=4 (766 aa)',
    'braf_historical_refseq': {'header': h1, 'length': len(s1),
                               'diff_vs_canonical_positions_1based': [i+1 for i in d1]},
    'braf_1992_mrna': {'header': h2, 'translated_length': len(t2),
                       'diff_vs_canonical_positions_1based': [i+1 for i in d2]},
    'note': ('The historical BRAF numbering (V599E, Davies 2002 Nature 417:949) '
             'shifts +1 to the canonical V600E because the historical reference '
             'lacks exactly one codon (3 nt) relative to the current canonical '
             'sequence. The difference and its position are recorded above; the '
             'alias is per-record evidence and is never generalized to other proteins.'),
}
(HERE / 'BRAF_HISTORICAL_EVIDENCE.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps(out, indent=1))
