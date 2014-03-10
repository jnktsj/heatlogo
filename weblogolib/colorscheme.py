
#  Copyright (c) 2003-2005 The Regents of the University of California.
#  Copyright (c) 2005 Gavin E. Crooks

#  This software is distributed under the MIT Open Source License.
#  <http://www.opensource.org/licenses/mit-license.html>
#
#  Permission is hereby granted, free of charge, to any person obtaining a 
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
#  THE SOFTWARE.

""" Popular color codings for nucleic and amino acids. 


Classes:
    ColorScheme -- A color scheme
    ColorGroup

    HeatScheme -- A color scheme for heatmap
    HeatGroup

o Color schemes for heatmap
    thermography
    heat_colors
    red_green
    red_blue
    magenta_cyan

o Color schemes without heatmap
    Generic
        monochrome

    Nucleotides
        nucleotide
        base pairing

    Amino Acid
        hydrophobicity
        chemistry
        charge
        taylor

Status : Beta - Needs documentation.

"""
# Good online references include bioruby and the JalView alignment editor.
# Clamp, M., Cuff, J., Searle, S. M. and Barton, G. J. (2004), 
# "The Jalview Java Alignment Editor," Bioinformatics, 12, 426-7
# http://www.jalview.org


from corebio import seq
from color import Color

class ColorScheme(object):
    """ A coloring of an alphabet.
    
    title : string            -- A human readable description
    defualt_color : Color           --
    groups : list of color groups 
    alphabet : string               -- The set of colored symbols
    color -- A map between a symbol and a Coloring
    

    """
    
    def __init__(self, 
                groups = [], 
                title = "", 
                description = "",
                default_color = "black", 
                alphabet = seq.generic_alphabet) :
        """  """
        self.title= title
        self.description = description
        self.default_color = Color.from_string(default_color)
        self.groups = groups
        self.alphabet = alphabet
            
        color = {}

        if alphabet in [seq.codon_dna_alphabet, seq.codon_rna_alphabet]:
            for cg in groups :
                color[cg.symbols] = cg.color
                if cg.symbols not in alphabet :
                    raise KeyError("Colored symbol does not exist in alphabet.")
        else:
            for cg in groups :
                for s in cg.symbols :
                    color[s] = cg.color
                    if s not in alphabet :
                        raise KeyError("Colored symbol does not exist in alphabet.")

        self._color = color

    def color(self, symbol) :
        if symbol in self._color :
            return self._color[symbol]
        return self.default_color


        
class ColorGroup(object) :
    """Associate a group of symbols with a color"""
    def __init__(self, symbols, color, description=None) :
        self.symbols = symbols              
        self.color =  Color.from_string(color)
        self.description = description

         
monochrome = ColorScheme([]) # This list intentionally left blank
               
# From makelogo
nucleotide = ColorScheme([
    ColorGroup("G", "orange"),
    ColorGroup("TU", "red"),
    ColorGroup("C",  "blue"),
    ColorGroup("A",  "green") 
    ]) 

base_pairing = ColorScheme([
    ColorGroup("TAU",  "darkorange", "Weak (2 Watson-Crick hydrogen bonds)"),
    ColorGroup("GC",    "blue", "Strong (3 Watson-Crick hydrogen bonds)")],
    )

# From Crooks2004c-Proteins-SeqStr.pdf
hydrophobicity = ColorScheme([
    ColorGroup( "RKDENQ",   "blue", "hydrophilic"),
    ColorGroup( "SGHTAP",   "green", "neutral"  ),
    ColorGroup( "YVMCLFIW", "black",  "hydrophobic") ],
    alphabet = seq.unambiguous_protein_alphabet
    )

# from makelogo
chemistry = ColorScheme([
  ColorGroup( "GSTYC",  "green",   "polar"),
  ColorGroup( "NQ",      "purple", "neutral"), 
  ColorGroup( "KRH",     "blue",   "basic"),
  ColorGroup( "DE",      "red",    "acidic"),
  ColorGroup("PAWFLIMV", "black",  "hydrophobic") ],
  alphabet = seq.unambiguous_protein_alphabet
  )   

charge = ColorScheme([
    ColorGroup("KRH", "blue", "Positive" ),
    ColorGroup( "DE", "red", "Negative") ],
    alphabet = seq.unambiguous_protein_alphabet
    )


taylor = ColorScheme([
    ColorGroup( 'A', '#CCFF00' ),
    ColorGroup( 'C', '#FFFF00' ),
    ColorGroup( 'D', '#FF0000'),
    ColorGroup( 'E', '#FF0066' ),
    ColorGroup( 'F', '#00FF66'),
    ColorGroup( 'G', '#FF9900'),
    ColorGroup( 'H', '#0066FF'),
    ColorGroup( 'I', '#66FF00'),
    ColorGroup( 'K', '#6600FF'),
    ColorGroup( 'L', '#33FF00'),
    ColorGroup( 'M', '#00FF00'),
    ColorGroup( 'N', '#CC00FF'),
    ColorGroup( 'P', '#FFCC00'),
    ColorGroup( 'Q', '#FF00CC'),
    ColorGroup( 'R', '#0000FF'),
    ColorGroup( 'S', '#FF3300'),
    ColorGroup( 'T', '#FF6600'),
    ColorGroup( 'V', '#99FF00'),
    ColorGroup( 'W', '#00CCFF'),
    ColorGroup( 'Y', '#00FFCC')],
    title = "Taylor",
    description = "W. Taylor, Protein Engineering, Vol 10 , 743-746 (1997)",
    alphabet = seq.unambiguous_protein_alphabet
    )


class HeatScheme(object):
    """ A coloring of a p-value.
    
    title : string
    groups : list of color groups 
    pvalue : float -- The set of colored pvalues
    color -- A map between a p-value and a coloring
    """    
    def __init__(self, groups = [], title = "") :
        self.title= title
        self.groups = groups
        color = {}
        for cg in groups:
            color[cg.pvalues] = cg.color
            # if float(cg.pvalues) < -1 or 1 < float(cg.pvalues):
            #     raise KeyError("P-value range should be in -1 < p < 1.")
        self._color = color

    def color(self, pvalue) :
        if pvalue in self._color :
            return self._color[pvalue]
        else:
            raise KeyError(str(pvalue)+" is not defined in heat color scheme.")

class HeatGroup(object) :
    """Associate a p-vaulue with a color"""
    def __init__(self, pvalues, color) :
        self.pvalues = pvalues
        self.color =  Color.from_string(color)

thermography = HeatScheme( [
    HeatGroup( '+0.01' , '#FF0000' ),
    HeatGroup( '+0.025', '#FF6600' ),
    HeatGroup( '+0.05' , '#FFCC00' ),
    HeatGroup( '+0.1'  , '#FFFF00' ),
    HeatGroup( '+0.3'  , '#CCFF33' ),
    HeatGroup( ' 0.5'  , '#00FF00' ),
    HeatGroup( '-0.3'  , '#00FFCC' ),
    HeatGroup( '-0.1'  , '#00FFFF' ),
    HeatGroup( '-0.05' , '#0066FF' ),
    HeatGroup( '-0.025', '#0000FF' ),
    HeatGroup( '-0.01' , '#330099' ) ],
    title = "Thermographic colors" )

red_green = HeatScheme( [
    HeatGroup( '+0.01' , '#FF0000' ),
    HeatGroup( '+0.025', '#FF2020' ),
    HeatGroup( '+0.05' , '#FF4040' ),
    HeatGroup( '+0.1'  , '#FF8080' ),
    HeatGroup( '+0.3'  , '#FFC0C0' ),
    HeatGroup( ' 0.5'  , '#E0FFFF' ),
    HeatGroup( '-0.3'  , '#C0FFC0' ),
    HeatGroup( '-0.1'  , '#80FF80' ),
    HeatGroup( '-0.05' , '#40FF40' ),
    HeatGroup( '-0.025', '#20FF20' ),
    HeatGroup( '-0.01' , '#00FF00' ) ],
    title = "Gradational colors between red and green" )

red_blue = HeatScheme( [
    HeatGroup( '+0.01' , '#FF0000' ),
    HeatGroup( '+0.025', '#FF2020' ),
    HeatGroup( '+0.05' , '#FF4040' ),
    HeatGroup( '+0.1'  , '#FF8080' ),
    HeatGroup( '+0.3'  , '#FFC0C0' ),
    HeatGroup( ' 0.5'  , '#E0FFFF' ),
    HeatGroup( '-0.3'  , '#C0C0FF' ),
    HeatGroup( '-0.1'  , '#8080FF' ),
    HeatGroup( '-0.05' , '#4040FF' ),
    HeatGroup( '-0.025', '#2020FF' ),
    HeatGroup( '-0.01' , '#0000FF' ) ],
    title = "Gradational colors between red and blue" )

magenta_cyan = HeatScheme( [
    HeatGroup( '+0.01' , '#FF00FF' ),
    HeatGroup( '+0.025', '#FF40FF' ),
    HeatGroup( '+0.05' , '#FF60FF' ),
    HeatGroup( '+0.1'  , '#FFA0FF' ),
    HeatGroup( '+0.3'  , '#FFE0FF' ),
    HeatGroup( ' 0.5'  , '#FFFFE0' ),
    HeatGroup( '-0.3'  , '#E0FFFF' ),
    HeatGroup( '-0.1'  , '#A0FFFF' ),
    HeatGroup( '-0.05' , '#60FFFF' ),
    HeatGroup( '-0.025', '#40FFFF' ),
    HeatGroup( '-0.01' , '#00FFFF' ) ],
    title = "Gradational colors between magenta and cyan" )


codons_dna = ColorScheme([
    ColorGroup( 'CAT', '#00FFFF'),
    ColorGroup( 'CAC', '#00FFFF'),
    ColorGroup( 'AAA', '#00FFFF'),
    ColorGroup( 'AAG', '#00FFFF'),
    ColorGroup( 'CGT', '#00FFFF'),
    ColorGroup( 'CGC', '#00FFFF'),
    ColorGroup( 'CGA', '#00FFFF'),
    ColorGroup( 'CGG', '#00FFFF'),
    ColorGroup( 'AGA', '#00FFFF'),
    ColorGroup( 'AGG', '#00FFFF'),
    ColorGroup( 'GAT', '#FF0000'),
    ColorGroup( 'GAC', '#FF0000'),
    ColorGroup( 'GAA', '#FF0000'),
    ColorGroup( 'GAG', '#FF0000'),
    ColorGroup( 'TCT', '#00FF00'),
    ColorGroup( 'TCC', '#00FF00'),
    ColorGroup( 'TCA', '#00FF00'),
    ColorGroup( 'TCG', '#00FF00'),
    ColorGroup( 'AGT', '#00FF00'),
    ColorGroup( 'AGC', '#00FF00'),
    ColorGroup( 'ACT', '#00FF00'),
    ColorGroup( 'ACC', '#00FF00'),
    ColorGroup( 'ACA', '#00FF00'),
    ColorGroup( 'ACG', '#00FF00'),
    ColorGroup( 'CAA', '#00FF00'),
    ColorGroup( 'CAG', '#00FF00'),
    ColorGroup( 'AAT', '#00FF00'),
    ColorGroup( 'AAC', '#00FF00'),
    ColorGroup( 'GCT', '#5555FF'),
    ColorGroup( 'GCC', '#5555FF'),
    ColorGroup( 'GCA', '#5555FF'),
    ColorGroup( 'GCG', '#5555FF'),
    ColorGroup( 'GTT', '#5555FF'),
    ColorGroup( 'GTC', '#5555FF'),
    ColorGroup( 'GTA', '#5555FF'),
    ColorGroup( 'GTG', '#5555FF'),
    ColorGroup( 'CTT', '#5555FF'),
    ColorGroup( 'CTC', '#5555FF'),
    ColorGroup( 'CTA', '#5555FF'),
    ColorGroup( 'CTG', '#5555FF'),
    ColorGroup( 'TTA', '#5555FF'),
    ColorGroup( 'TTG', '#5555FF'),
    ColorGroup( 'ATT', '#5555FF'),
    ColorGroup( 'ATC', '#5555FF'),
    ColorGroup( 'ATA', '#5555FF'),
    ColorGroup( 'ATG', '#5555FF'),
    ColorGroup( 'TTT', '#FF00FF'),
    ColorGroup( 'TTC', '#FF00FF'),
    ColorGroup( 'TAT', '#FF00FF'),
    ColorGroup( 'TAC', '#FF00FF'),
    ColorGroup( 'TGG', '#FF00FF'),
    ColorGroup( 'GGT', '#996600'),
    ColorGroup( 'GGC', '#996600'),
    ColorGroup( 'GGA', '#996600'),
    ColorGroup( 'GGG', '#996600'),
    ColorGroup( 'CCT', '#996600'),
    ColorGroup( 'CCC', '#996600'),
    ColorGroup( 'CCA', '#996600'),
    ColorGroup( 'CCG', '#996600'),
    ColorGroup( 'TGT', '#FFFF00'),
    ColorGroup( 'TGC', '#FFFF00'),
    ColorGroup( 'TAA', '#000000'),
    ColorGroup( 'TAG', '#000000'),
    ColorGroup( 'TGA', '#000000')],  
    alphabet = seq.codon_dna_alphabet
    )

codons_rna = ColorScheme([
    ColorGroup( 'CAU', '#00FFFF'),
    ColorGroup( 'CAC', '#00FFFF'),
    ColorGroup( 'AAA', '#00FFFF'),
    ColorGroup( 'AAG', '#00FFFF'),
    ColorGroup( 'CGU', '#00FFFF'),
    ColorGroup( 'CGC', '#00FFFF'),
    ColorGroup( 'CGA', '#00FFFF'),
    ColorGroup( 'CGG', '#00FFFF'),
    ColorGroup( 'AGA', '#00FFFF'),
    ColorGroup( 'AGG', '#00FFFF'),
    ColorGroup( 'GAU', '#FF0000'),
    ColorGroup( 'GAC', '#FF0000'),
    ColorGroup( 'GAA', '#FF0000'),
    ColorGroup( 'GAG', '#FF0000'),
    ColorGroup( 'UCU', '#00FF00'),
    ColorGroup( 'UCC', '#00FF00'),
    ColorGroup( 'UCA', '#00FF00'),
    ColorGroup( 'UCG', '#00FF00'),
    ColorGroup( 'AGU', '#00FF00'),
    ColorGroup( 'AGC', '#00FF00'),
    ColorGroup( 'ACU', '#00FF00'),
    ColorGroup( 'ACC', '#00FF00'),
    ColorGroup( 'ACA', '#00FF00'),
    ColorGroup( 'ACG', '#00FF00'),
    ColorGroup( 'CAA', '#00FF00'),
    ColorGroup( 'CAG', '#00FF00'),
    ColorGroup( 'AAU', '#00FF00'),
    ColorGroup( 'AAC', '#00FF00'),
    ColorGroup( 'GCU', '#5555FF'),
    ColorGroup( 'GCC', '#5555FF'),
    ColorGroup( 'GCA', '#5555FF'),
    ColorGroup( 'GCG', '#5555FF'),
    ColorGroup( 'GUU', '#5555FF'),
    ColorGroup( 'GUC', '#5555FF'),
    ColorGroup( 'GUA', '#5555FF'),
    ColorGroup( 'GUG', '#5555FF'),
    ColorGroup( 'CUU', '#5555FF'),
    ColorGroup( 'CUC', '#5555FF'),
    ColorGroup( 'CUA', '#5555FF'),
    ColorGroup( 'CUG', '#5555FF'),
    ColorGroup( 'UUA', '#5555FF'),
    ColorGroup( 'UUG', '#5555FF'),
    ColorGroup( 'AUU', '#5555FF'),
    ColorGroup( 'AUC', '#5555FF'),
    ColorGroup( 'AUA', '#5555FF'),
    ColorGroup( 'AUG', '#5555FF'),
    ColorGroup( 'UUU', '#FF00FF'),
    ColorGroup( 'UUC', '#FF00FF'),
    ColorGroup( 'UAU', '#FF00FF'),
    ColorGroup( 'UAC', '#FF00FF'),
    ColorGroup( 'UGG', '#FF00FF'),
    ColorGroup( 'GGU', '#996600'),
    ColorGroup( 'GGC', '#996600'),
    ColorGroup( 'GGA', '#996600'),
    ColorGroup( 'GGG', '#996600'),
    ColorGroup( 'CCU', '#996600'),
    ColorGroup( 'CCC', '#996600'),
    ColorGroup( 'CCA', '#996600'),
    ColorGroup( 'CCG', '#996600'),
    ColorGroup( 'UGU', '#FFFF00'),
    ColorGroup( 'UGC', '#FFFF00'),
    ColorGroup( 'UAA', '#000000'),
    ColorGroup( 'UAG', '#000000'),
    ColorGroup( 'UGA', '#000000')],
    alphabet = seq.codon_rna_alphabet
    )
