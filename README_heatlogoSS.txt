HeatLogoSS
==========

HeatLogoSS generates sequence logos from biological sequence alignments
and visualizes sequence conservation with color-coded symbols representing
over-/under-representation of symbols (HeatLogo), and protein structure
information (i.e. secondary structure and disorder prediction; StructureLogo).


Requirement
-----------

MultiLogo requires the following environment to run:

  1) Python 2.x (x >= 5)
  2) Python array package 'numpy':
  3) Python scientific package 'scipy':             http://www.scipy.org/Download
  4) Secondary structure prediction tool 'PsiPred': http://bioinf.cs.ucl.ac.uk/psipred
  5) Disorder prediction tool 'Poodle-L':           http://mbs.cbrc.jp/poodle/poodle-l.html

  * Please replace 'runpsipred' source code in PsiPred package with the one
    under 'sslogolib' directory.  After replacing the source code, please
    set pathes of database, software etc. by following the PsiPred instruction.


Usage and Examples
------------------

For help on the command line interface, run:
    $ ./heatlogoSS --help
or simply type:
    $ ./heatlogoSS -h

To build a simple logo, execute:
    $ ./heatlogoSS --poodlel PATH --psipred PATH < aa-seqs.fa

'--poodlel' and '--psipred' are required options to run heatlogoSS.
For obtaining more promising prediction results, using whole protein
alignments as an input is recommended.

In addition to this, if 'montage' is installed in the computer,
heatlogoSS will merge HeatLogo and StructureLogo into one file.

If your sequences are very long and want to extract a specific
regions, run:
    $ ./heatlogoSS --poodlel PATH --psipred PATH \
                   --start 5 --end 15 < aa-seqs.fa

Like HeatLogo, HeatLogoSS can also generate position weight matrices of
P-values and probabilities of input sequences and their protein structures.
    $ ./heatlogoSS --poodlel PATH --psipred PATH \
                   --aa-pwm-pval FILENAME --aa-pwm-prob FILENAME \
                   --st-pwm-prob FILENAME < aa-seqs.fa


heatlogoSS has several common options between HeatLogo and StrucutureLogo.

--[*]-title      : Title of sequence and structure logo
--[*]-show-xaxis : Display X-axis in sequence and structure logo
--[*]-show-yaxis : Display Y-axis in sequence and structure logo
--[*]-yaxis      : Hight of Y-axis in units
--[*]-xlabel     : X-axis label
--[*]-ylabel     : Y-axis label
--[*]-ticmarks   : Distance between ticmarks

[*] can be replaced with 'aa' or 'st' for generating HeatLogo or StructureLogo.


For more detail, please type: ./heatlogoSS --help
