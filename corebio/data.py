#  Copyright (c) 2006, The Regents of the University of California, through 
#  Lawrence Berkeley National Laboratory (subject to receipt of any required
#  approvals from the U.S. Dept. of Energy).  All rights reserved.

#  This software is distributed under the new BSD Open Source License.
#  <http://www.opensource.org/licenses/bsd-license.html>
#
#  Redistribution and use in source and binary forms, with or without 
#  modification, are permitted provided that the following conditions are met: 
#
#  (1) Redistributions of source code must retain the above copyright notice, 
#  this list of conditions and the following disclaimer. 
#
#  (2) Redistributions in binary form must reproduce the above copyright 
#  notice, this list of conditions and the following disclaimer in the 
#  documentation and or other materials provided with the distribution. 
#
#  (3) Neither the name of the University of California, Lawrence Berkeley 
#  National Laboratory, U.S. Dept. of Energy nor the names of its contributors 
#  may be used to endorse or promote products derived from this software 
#  without specific prior written permission. 
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
#  POSSIBILITY OF SUCH DAMAGE. 

"""
Standard data used in computational biology.


To convert a property dictionary to a list :
>>> comp = [ amino_acid_composition[k] for k in amino_acid_letters]



Resources: 
    Various standard data files are included in the corebio distribution. These
    may be loaded with the data_string, data_stream or data_filename methods.
    A complete set of names is stored in 'resource_names'
 
BLOSUM Scoring Matrices
    Source: ftp://ftp.ncbi.nih.gov/repository/blocks/unix/blosum
    These are all new blast style with 1/3 bit scaling
    - blosum35
    - blosum45    
    - blosum62    
    - blosum40    
    - blosum50    
    - blosum80    
    - blosum100   

Other substitution scoring matrices:
    - dist20_comp 
    - pam250
    - pam120
    - vtml160
    
Description of database cross references :
    - dbxref.txt (http://www.expasy.org/cgi-bin/lists?dbxref.txt)

    
Attributes:
    - amino_acid_letters
        -- Standard codes for the 20 canonical amino acids, in alphabetic
        order.
        
    - amino_acid_alternative_letters
        -- Amino acid one letter codes, alphabetic by three letter codes.

    - amino_acid_extended_letters

    - dna_letters

    - dna_extended_letters

    - rna_letters
    
    - rna_extended_letters

    - dna_ambiguity 

    - rna_ambiguity
    
    - amino_acid_ambiguity
    
    - amino_acid_mass
        -- Monomer isotopically averaged molecular mass 
    
    - dna_mass
    
    - rna_mass
        
    - one_to_three      
        -- Map from standard 1 letter amino acid codes to standard three
        letter codes. 
        Ref: http://www.ebi.ac.uk/RESID/faq.html
      
    - standard_three_to_one
        -- Map from standard 3 letter amino acid codes to standard 1
        letter codes.
         
    - extended_three_to_one
        -- Map between three letter amino acid codes (first letter capitalized) 
        and standard one letter codes. This map contains many nonstandard three
        letter codes, used, for example, to specify chemically modified amino
        acids in PDB files.
        Ref: http://astral.berkeley.edu/ 
        Ref: http://www.ebi.ac.uk/RESID/faq.html

    - amino_acid_names

    - amino_acid_composition
        -- Average amino acid composition of proteins.
        Ref: McCaldon P., Argos P. Proteins 4:99-122 (1988).

    - kyte_doolittle_hydrophobicity 
        -- Kyte-Doolittle hydrophobicity scale.
        Ref: Kyte J., Doolittle R.F. J. Mol. Biol. 157:105-132 (1982)
        
    - nucleotide_names
    
    - amino_acid_accesible_surface_area
        -- Nominal maximum solvent accessoble area for unmodified amino acids,
        in square Angstroms.
        Ref: Sander & Rost, (1994), Proteins, 20:216-226


Status: Beta (Data needs to be proof checked.)    
"""

# FIXME: Proof check data
# FIXME: Add __all__

# The ExPasy ProtScale tool is a great source of amino acid properties.
# http://au.expasy.org/cgi-bin/protscale.pl       

from StringIO import StringIO
from corebio.utils import resource_string, resource_stream,resource_filename
import utils

# Explicitly list set of available data resources. We want to be able to access
# these resources in, for example, a webapp, without inadvertently allowing
# unrestricted read access to the local file system.

resource_names = [
    'blosum35',
    'blosum45',    
    'blosum62',    
    'blosum40',    
    'blosum50',    
    'blosum80',    
    'blosum100',   
    'dist20_comp', 
    'pam250',
    'pam120', 
    'dbxref.txt',
    'vtml160',
    ]

_resource_filenames = {
    'blosum35':    'data/blosum35.mat',
    'blosum45':    'data/blosum45.mat',    
    'blosum62':    'data/blosum62.mat',    
    'blosum40':    'data/blosum40.mat',    
    'blosum50':    'data/blosum50.mat',    
    'blosum80':    'data/blosum80.mat',    
    'blosum100':   'data/blosum100.mat',   
    'dist20_comp': 'data/dist20_comp.mat', 
    'pam250':      'data/pam250.mat',
    'pam120':      'data/pam120.mat', 
    'dbxref.txt' : 'data/dbxref.txt',
    'vtml160' :    'data/vtml160',
    }

# TODO: Substitution matrix parser, SeqMatrix.read
# _resource_parsers = {}

def data_string( name ): 
    """Load the specified resource as a string."""
    fn = _resource_filenames[name]
    return resource_string(__name__, fn , __file__)    

def data_stream( name ):
    """Provide an open file handle to the specified resource."""
    fn = _resource_filenames[name]
    return resource_stream(__name__, fn , __file__)    

def data_filename( name ): 
    """Provide a filename for the given resource in the local filesystem."""
    fn = _resource_filenames[name]
    return resource_filename(__name__, fn, __file__)            

#def data_object( name, parser = None) :
#    if parser is None : 
#        if name in _resource_parsers :
#            parser = _resource_parsers[name]
#        else :
#            parser = str    
#    return parser( data_stream(name) )


amino_acid_letters = "ACDEFGHIKLMNPQRSTVWY"

amino_acid_alternative_letters = "ARNDCQEGHILKMFPSTWYV"

amino_acid_extended_letters = "ACDEFGHIKLMNOPQRSTUVWYBJZX*-"


dna_letters = "GATC"
dna_extended_letters = "GATCRYWSMKHBVDN"

rna_letters = "GAUC"  
rna_extended_letters = "GAUCRYWSMKHBVDN"


dna_ambiguity = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "X": "GATC",
    "N": "GATC",
}

rna_ambiguity = {
    "A": "A",
    "C": "C",
    "G": "G",
    "U": "U",
    "M": "AC",
    "R": "AG",
    "W": "AU",
    "S": "CG",
    "Y": "CU",
    "K": "GU",
    "V": "ACG",
    "H": "ACU",
    "D": "AGU",
    "B": "CGU",
    "X": "GAUC",
    "N": "GAUC",
}

amino_acid_ambiguity = {
    "A": "A",
    "B": "ND",
    "C": "C",
    "D": "D",
    "E": "E",
    "F": "F",
    "G": "G",
    "H": "H",
    "I": "I",
    "K": "K",
    "L": "L",
    "M": "M",
    "N": "N",
    "P": "P",
    "Q": "Q",
    "R": "R",
    "S": "S",
    "T": "T",
    "V": "V",
    "W": "W",
    "X": "ACDEFGHIKLMNPQRSTVWY",
    "Y": "Y",
    "Z": "QE",
    "J": "IL",
    'U': 'U',
    'O': 'O',
}


# Monomer isotopically averaged molecular mass 
# Data Checked GEC Nov 2006
amino_acid_mass = {
    "A": 89.09,
    "B" : 132.66,  # Averaged proportional to amino_acid_composition
    "C": 121.16,
    "D": 133.10,
    "E": 147.13,
    "F": 165.19,
    "G": 75.07,
    "H": 155.16, 
    "I": 131.18,
    "J": 131.18,
    "K": 146.19,
    "L": 131.18,
    "M": 149.21,
    "N": 132.12,
    # "O" : ???, # TODO
    "P": 115.13,
    "Q": 146.15,
    "R": 174.20,
    "S": 105.09,
    "T": 119.12,
    "U" : 168.05,
    "V": 117.15,
    "W": 204.23,
    "X" : 129.15, # Averaged proportional to amino_acid_composition  
    "Y": 181.19,
    "Z" : 146.76, # Averaged proportional to amino_acid_composition    
    }
    
dna_mass = {
    "A": 347.,
    "C": 323.,
    "G": 363.,
    "T": 322.,
    }

rna_mass = {
    "A": 363.,
    "C": 319.,
    "G": 379.,
    "U": 340.,
}

one_to_three = {
    'A':'Ala', 'B':'Asx', 'C':'Cys', 'D':'Asp',
    'E':'Glu', 'F':'Phe', 'G':'Gly', 'H':'His',
    'I':'Ile', 'K':'Lys', 'L':'Leu', 'M':'Met',
    'N':'Asn', 'P':'Pro', 'Q':'Gln', 'R':'Arg',
    'S':'Ser', 'T':'Thr', 'V':'Val', 'W':'Trp',
    'Y':'Tyr', 'Z':'Glx', 'X':'Xaa', 
    'U':'Sec', 'J':'Xle', 'O':'Pyl'
    }


standard_three_to_one = utils.invert_dict(one_to_three)

extended_three_to_one= {
'2as':'D', '3ah':'H', '5hp':'E', 'Acl':'R', 'Agm':'R', 'Aib':'A', 'Ala':'A', 'Alm':'A', 'Alo':'T', 'Aly':'K', 'Arg':'R', 'Arm':'R', 'Asa':'D', 'Asb':'D', 'Ask':'D', 'Asl':'D', 'Asn':'N', 'Asp':'D', 'Asq':'D', 'Asx':'B', 'Aya':'A', 'Bcs':'C', 'Bhd':'D', 'Bmt':'T', 'Bnn':'A', 'Buc':'C', 'Bug':'L', 'C5c':'C', 'C6c':'C', 'Ccs':'C', 'Cea':'C', 'Cgu':'E', 'Chg':'A', 'Cle':'L', 'Cme':'C', 'Csd':'A', 'Cso':'C', 'Csp':'C', 'Css':'C', 'Csw':'C', 'Csx':'C', 'Cxm':'M', 'Cy1':'C', 'Cy3':'C', 'Cyg':'C', 'Cym':'C', 'Cyq':'C', 'Cys':'C', 'Dah':'F', 'Dal':'A', 'Dar':'R', 'Das':'D', 'Dcy':'C', 'Dgl':'E', 'Dgn':'Q', 'Dha':'A', 'Dhi':'H', 'Dil':'I', 'Div':'V', 'Dle':'L', 'Dly':'K', 'Dnp':'A', 'Dpn':'F', 'Dpr':'P', 'Dsn':'S', 'Dsp':'D', 'Dth':'T', 'Dtr':'W', 'Dty':'Y', 'Dva':'V', 'Efc':'C', 'Fla':'A', 'Fme':'M', 'Ggl':'E', 'Gl3':'G', 'Gln':'Q', 'Glu':'E', 'Glx':'Z', 'Gly':'G', 'Glz':'G', 'Gma':'E', 'Gsc':'G', 'Hac':'A', 'Har':'R', 'Hic':'H', 'Hip':'H', 'His':'H', 'Hmr':'R', 'Hpq':'F', 'Htr':'W', 'Hyp':'P', 'Iil':'I', 'Ile':'I', 'Iyr':'Y', 'Kcx':'K', 'Leu':'L', 'Llp':'K', 'Lly':'K', 'Ltr':'W', 'Lym':'K', 'Lys':'K', 'Lyz':'K', 'Maa':'A', 'Men':'N', 'Met':'M', 'Mhs':'H', 'Mis':'S', 'Mle':'L', 'Mpq':'G', 'Msa':'G', 'Mse':'M', 'Mva':'V', 'Nem':'H', 'Nep':'H', 'Nle':'L', 'Nln':'L', 'Nlp':'L', 'Nmc':'G', 'Oas':'S', 'Ocs':'C', 'Omt':'M', 'Paq':'Y', 'Pca':'E', 'Pec':'C', 'Phe':'F', 'Phi':'F', 'Phl':'F', 'Pr3':'C', 'Pro':'P', 'Prr':'A', 'Ptr':'Y', 'Pyl':'O', 'Sac':'S', 'Sar':'G', 'Sch':'C', 'Scs':'C', 'Scy':'C', 'Sec':'U', 'Sel':'U', 'Sep':'S', 'Ser':'S', 'Set':'S', 'Shc':'C', 'Shr':'K', 'Smc':'C', 'Soc':'C', 'Sty':'Y', 'Sva':'S', 'Ter':'*', 'Thr':'T', 'Tih':'A', 'Tpl':'W', 'Tpo':'T', 'Tpq':'A', 'Trg':'K', 'Tro':'W', 'Trp':'W', 'Tyb':'Y', 'Tyq':'Y', 'Tyr':'Y', 'Tys':'Y', 'Tyy':'Y', 'Unk':'X', 'Val':'V', 'Xaa':'X', 'Xer':'X', 'Xle':'J'}
# Initial table is from the ASTRAL RAF release notes.
# added UNK
# Extra IUPAC: Xle, Xaa, Sec, Pyl
# The following have been seen in biopython code.
# Ter : '*'     Termination
# Sel : 'U'     A typo for Sec, selenocysteine? 
# Xer : 'X'     Another alternative for unknown?


amino_acid_names = {
    'A'	: 'alanine',	
    'M'	: 'methionine',  
    'C'	: 'cysteine',
    'N'	: 'asparagine',
    'D'	: 'aspartic acid',
    'P'	: 'proline',
    'E'	: 'glutamic acid',
    'Q'	: 'glutamine',
    'F'	: 'phenylalanine',
    'R'	: 'arginine',
    'G'	: 'glycine',	
    'S'	: 'serine',
    'H'	: 'histidine',	
    'T' : 'threonine',
    'I'	: 'isoleucine',	
    'V'	: 'valine',
    'K'	: 'lysine',
    'W'	: 'tryptophan', 
    'L'	: 'leucine',	
    'Y'	: 'tyrosine', 
    'B' : 'aspartic acid or asparagine',
    'J' : 'leucine or isoleucine',
    'X' : 'unknown',
    'Z' : 'glutamic acid or glutamine',
    'U' : 'selenocysteine',
    'O' : 'pyrrolysine',
    '*' : 'translation stop',
    '-' : 'gap'
    }

amino_acid_composition = dict(
    A = .082, R = .057, N = .044, D = .053, C = .017, 
    Q = .040, E = .062, G = .072, H = .022, I = .052,  
    L = .090, K = .057, M = .024, F =.039, P = .051, 
    S = .069, T = .058, W = .013, Y= .032, V =.066 )


kyte_doolittle_hydrophobicity = dict(
    A=1.8, R=-4.5, N=-3.5, D=-3.5,  C=2.5, 
    Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
    L=3.8, K=-3.9,  M=1.9,  F=2.8, P=-1.6,
    S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2 )


nucleotide_names = { 
    'A' : 'Adenosine',
    'C'	: 'Cytidine',
    'G'	: 'Guanine',
    'T'	: 'Thymidine',
    'U'	: 'Uracil',
    'R'	: 'G A (puRine)',
    'Y'	: 'T C (pYrimidine)',
    'K'	: 'G T (Ketone)',
    'M'	: 'A C (aMino group)',
    'S'	: 'G C (Strong interaction)',
    'W'	: 'A T (Weak interaction)',
    'B'	: 'G T C (not A) (B comes after A)',
    'D'	: 'G A T (not C) (D comes after C)',
    'H'	: 'A C T (not G) (H comes after G)',
    'V'	: 'G C A (not T, not U) (V comes after U)',
    'N' : 'A G C T (aNy)',
    '-' : 'gap', 
    }
    

# TODO: CHECK VALUES, UNITS    
amino_acid_accesible_surface_area = {
    'A' : 106.0,
    'C' : 135.0,
    'D' : 163.0,
    'E' : 194.0,
    'F' : 197.0,
    'G' : 84.0,
    'H' : 184.0,
    'I' : 169.0,
    'K' : 205.0,
    'L' : 164.0,
    'M' : 188.0,
    'N' : 157.0,
    'P' : 136.0,
    'Q' : 198.0,
    'R' : 248.0,
    'S' : 130.0,
    'T' : 142.0,
    'V' : 142.0,
    'W' : 227.0,
    'Y' : 222.0
    }


ecoli_codon_composition = dict(
    CTT = 0.7616, ATG = 1.5872, ACA = 0.4096, ACG = 0.7360,
    ATC = 1.1648, AAC = 1.5616, ATA = 0.2368, AGG = 0.1024,
    CCT = 0.5376, ACT = 0.5120, AGC = 1.0624, AAG = 0.7744,
    AGA = 0.0896, CAT = 1.0112, AAT = 1.4016, ATT = 1.9520,
    CTG = 3.0016, CTA = 0.3392, CTC = 0.6720, CAC = 0.8384,
    AAA = 2.1248, CCG = 1.7088, AGT = 0.4608, CCA = 0.4224,
    CAA = 0.7744, CCC = 0.4096, TAT = 1.0752, GGT = 1.3632,
    TGT = 0.3776, CGA = 0.2752, CAG = 1.7728, TCT = 0.3648,
    GAT = 2.4256, CGG = 0.2624, TTT = 1.2608, TGC = 0.5120,
    GGG = 0.5504, TAG = 1e-06 , GGA = 0.5888, TAA = 0.1152,
    GGC = 2.1376, TAC = 0.9344, TTC = 0.9600, TCG = 0.5120,
    TTA = 0.9728, TTG = 0.7616, TCC = 0.3520, ACC = 1.4592,
    TCA = 0.4992, GCA = 1.3504, GTA = 0.7360, GCC = 2.0224,
    GTC = 0.7488, GCG = 2.4640, GTG = 1.6896, GAG = 1.1776,
    GTT = 1.0752, GCT = 0.6848, TGA = 0.0640, GAC = 1.3120,
    CGT = 1.3504, TGG = 0.6848, GAA = 2.7968, CGC = 1.664 )

human_codon_composition = dict(
    CTT = 0.8448, ATG = 1.4080, ACA = 0.9664, ACG = 0.3904,
    ATC = 1.3312, AAC = 1.2224, ATA = 0.4800, AGG = 0.7680,
    CCT = 1.1200, ACT = 0.8384, AGC = 1.2480, AAG = 2.0416,
    AGA = 0.7808, CAT = 0.6976, AAT = 1.0880, ATT = 1.0240,
    CTG = 2.5344, CTA = 0.4608, CTC = 1.2544, CAC = 0.9664,
    AAA = 1.5616, CCG = 0.4416, AGT = 0.7744, CCA = 1.0816,
    CAA = 0.7872, CCC = 1.2672, TAT = 0.7808, GGT = 0.6912,
    TGT = 0.6784, CGA = 0.3968, CAG = 2.1888, TCT = 0.9728,
    GAT = 1.3952, CGG = 0.7296, TTT = 1.1264, TGC = 0.8064,
    GGG = 1.0560, TAG = 0.0512, GGA = 1.0560, TAA = 0.0640,
    GGC = 1.4208, TAC = 0.9792, TTC = 1.2992, TCG = 0.2816,
    TTA = 0.4928, TTG = 0.8256, TCC = 1.1328, ACC = 1.2096,
    TCA = 0.7808, GCA = 1.0112, GTA = 0.4544, GCC = 1.7728,
    GTC = 0.9280, GCG = 0.4736, GTG = 1.7984, GAG = 2.5344,
    GTT = 0.7040, GCT = 1.1776, TGA = 0.1024, GAC = 1.6064,
    CGT = 0.2880, TGG = 0.8448, GAA = 1.8560, CGC = 0.6656 )

yeast_codon_composition = dict(
    CTT = 0.7872, ATG = 1.3376, ACA = 1.1392, ACG = 0.5120,
    ATC = 1.1008, AAC = 1.5872, ATA = 1.1392, AGG = 0.5888,
    CCT = 0.8640, ACT = 1.2992, AGC = 0.6272, AAG = 1.9712,
    AGA = 1.3632, CAT = 0.8704, AAT = 2.2848, ATT = 1.9264,
    CTG = 0.6720, CTA = 0.8576, CTC = 0.3456, CAC = 0.4992,
    AAA = 2.6816, CCG = 0.3392, AGT = 0.9088, CCA = 1.1712,
    CAA = 1.7472, CCC = 0.4352, TAT = 1.2032, GGT = 1.5296,
    TGT = 0.5184, CGA = 0.1920, CAG = 0.7744, TCT = 1.5040,
    GAT = 2.4064, CGG = 0.1088, TTT = 1.6704, TGC = 0.3072,
    GGG = 0.3840, TAG = 0.0320, GGA = 0.6976, TAA = 0.0704,
    GGC = 0.6272, TAC = 0.9472, TTC = 1.1776, TCG = 0.5504,
    TTA = 1.6768, TTG = 1.7408, TCC = 0.9088, ACC = 0.8128,
    TCA = 1.1968, GCA = 1.0368, GTA = 0.7552, GCC = 0.8064,
    GTC = 0.7552, GCG = 0.3968, GTG = 0.6912, GAG = 1.2288,
    GTT = 1.4144, GCT = 1.3568, TGA = 0.0448, GAC = 1.2928,
    CGT = 0.4096, TGG = 0.6656, GAA = 2.9184, CGC = 0.1664 )
