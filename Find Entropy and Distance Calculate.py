import math
import collections

def distance(a, b):
    if a == b:
        return 0
    elif a == "-" or b == "-":
        return 1
    else:
        return 2

def w_sum(a, b, w=1):
    if a == b:
        return 3 * w
    elif a == "-" or b == "-":
        return -1 * w
    else:
        return 2 * w

seq1 = 'A-TTG-CTT'
seq2 = 'ACCTGGC-T'
seq3 = 'ACTA-GCTA'
seq4 = 'AGTAGCATT'

list(zip(seq1, seq2, seq3, seq4)) == [('A', 'A', 'A', 'A'),
                                      ('-', 'C', 'C', 'G'),
                                      ('T', 'C', 'T', 'T'),
                                      ('T', 'T', 'A', 'A'),
                                      ('G', 'G', '-', 'G'),
                                      ('-', 'G', 'G', 'C'),
                                      ('C', 'C', 'C', 'A'),
                                      ('T', '-', 'T', 'T'),
                                      ('T', 'T', 'A', 'T')]

from itertools import combinations

list(combinations(('-', 'C', 'C', 'G'), 2)) == [('-', 'C'),
                                                ('-', 'C'),
                                                ('-', 'G'),
                                                ('C', 'C'),
                                                ('C', 'G'),
                                                ('C', 'G')]

consensus = 0
sop = 0
for position in zip(seq1, seq2, seq3, seq4):
    for pair in combinations(position, 2):
        consensus+= distance(*pair)
        sop += w_sum(*pair)

print("Of the gene sequence given in this section. The sum of pairs score and consensus score are calculated.")
print("Sums of Pair Score: ",sop)
print("Consensus Score: ",consensus)


def entropy(sequence):
    m=len(sequence)
    bases=collections.Counter([tmp_base for tmp_base in sequence])
    entropy_value=0
    for base in bases:
        n_i=bases[base]

        p_i=n_i / float(m)
        entropy_i=p_i * (math.log(p_i, 2))
        entropy_value+=entropy_i
    return entropy_value * -1

print("We find the entropy of the gene sequence given in this section.")

print("A-TTG-CTT Entropi: ",entropy("A-TTG-CTT"))
print("ACCTGGC-T Entropi: ",entropy("ACCTGGC-T"))
print("ACTA-GCTA Entropi: ",entropy("ACTA-GCTA"))
print("AGTAGCATT Entropi: ",entropy("AGTAGCATT"))
