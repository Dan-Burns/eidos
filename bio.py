import numpy as np

# from biopython register_ncbi_table (no stop codons)
codon_table={
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y",                          
        "TGT": "C", "TGC": "C",             "TGG": "W",   
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }



'''
Atchley, William R., Jieping Zhao, Andrew D. Fernandes, and Tanja Drüke. 2005. 
“Solving the Protein Sequence Metric Problem.” 
Proceedings of the National Academy of Sciences of the United States of America
 102 (18): 6395–6400.
'''
# five factor solution scores for 54 amino acid attributes
atchley_factors = {
    'A': np.array([-0.591, -1.302, -0.733, 1.570, -0.146]),
    'C': np.array([-1.343, 0.465, -0.862, -1.020, -0.255]),
    'D': np.array([1.050, 0.302, -3.656, -0.259, -3.242]),
    'E': np.array([1.357, -1.453, 1.477, 0.113, -0.837]),
    'F': np.array([-1.006, -0.590, 1.891, -0.397, 0.412]),
    'G': np.array([-0.384, 1.652, 1.330, 1.045, 2.064]),
    'H': np.array([0.336, -0.417, -1.673, -1.474, -0.078]),
    'I': np.array([-1.239, -0.547, 2.131, 0.393, 0.816]),
    'K': np.array([1.831, -0.561, 0.533, -0.277, 1.648]),
    'L': np.array([-1.019, -0.987, -1.505, 1.266, -0.912]),
    'M': np.array([-0.663, -1.524, 2.219, -1.005, 1.212]),
    'N': np.array([0.945, 0.828, 1.299, -0.169, 0.933]),
    'P': np.array([0.189, 2.081, -1.628, 0.421, -1.392]),
    'Q': np.array([0.931, -0.179, -3.005, -0.503, -1.853]),
    'R': np.array([1.538, -0.055, 1.502, 0.440, 2.897]),
    'S': np.array([-0.228, 1.399, -4.760, 0.670, -2.647]),
    'T': np.array([-0.032, 0.326, 2.213, 0.908, 1.313]),
    'V': np.array([-1.337, -0.279, -0.544, 1.242, -1.262]),
    'W': np.array([-0.595, 0.009, 0.672, -2.128, -0.184]),
    'Y': np.array([0.260, 0.830, 3.097, -0.838, 1.512])
   }

# names of the distal side chain atoms
distal_atoms = {'MET': 'CE',
               'LEU': 'CD1',
               'ALA': 'CB',
               'LYS': 'NZ',
               'ARG': 'NH1',
               'ILE': 'CD1',
               'CYS': 'SG',
               'ASP': 'OD1',
               'VAL': 'CG1',
               'GLY': 'HA2',
               'THR': 'OG1',
               'ASN': 'OD1',
               'PHE': 'CE1', # ortho position
               'GLU': 'OE1',
               'SER': 'OG',
               'PRO': 'CD',
               'TYR': 'CE1', # ortho position - else 'OH'
               'GLN': 'OE1',
               'HIS': 'NE2',
               'TRP': 'CH2',
               }