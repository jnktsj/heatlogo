========
HeatLogo
========
version: 1.0
.............

General description
-------------------

HeatLogo is a tool for generating sequence logos from biological sequence
alignments and for coloring each symbol in the logos according to its P-value.
It can be run on the command line, as a standalone webserver, as a CGI webapp,
or as a python library.

Requirement
-----------

HeatLogo requires the following environment to run:

- Python 2.x (x >= 5)
- Python array package 'numpy'
- Python scientific package 'scipy'

('numpy' and 'scipy' can be freely downloaded from: http://www.scipy.org/Download)

Usage
-----
Create sequence logos from biological sequence alignments and color symbols
according to P-values::

    $ heatlogo [options] < sequence_data.fa > sequence_logo.eps

Options
-------

  -h, --help
      Show all options, their default settings and exit

Input/Output Options
~~~~~~~~~~~~~~~~~~~~

  -i, --input FILENAME         Sequence input file (default: stdin)

  -f, --input-format FORMAT    Type of multiple sequence alignment or position weight matrix
                               file: (clustal, fasta, plain, msf, genbank, nbrf, nexus phylip
                               stockholm, intelligenetics, table, array, transfac)

  -o, --output FILENAME        Output file (default: stdout)

  -F, --output-format FORMAT   Format of output: eps (default), png, png_print, pdf, jpeg, svg, logodata.

      --pwm-prob FILENAME      Output position weight matrix composed of logarithmic probabilities

      --pwm-pval FILENAME      Output position weight matrix of P-values

Logo Data Options
~~~~~~~~~~~~~~~~~
  -A, --sequence-type TYPE     The type of sequence data: 'protein', 'rna', 'dna', or 'codon'
                               (See also: --codon-frame)

  -a, --alphabet ALPHABET      The set of symbols to count, e.g. 'AGTC'. All characters not in
                               the alphabet are ignored. If neither the alphabet nor sequence-type
                               are specified then heatlogo will examine the input data and make
                               an educated guess. (See also: --sequence-type, --ignore-lower-case)

  -R, --codon-frame NUMBER     Codon reading frame (default: +1): [+1, +2, +3] indicate the forward
                               strand frame, and [-1, -2, -3] indicate the reverse strand frame

  -b, --composition COMP       The expected composition of the sequences. HeatLogo accepts four 
                               types of formats:

                                   1) String: 'auto' (default), 'equiprobable', or 'none'
                                      (do not perform any compositional adjustment, Heatmap: off).
                                      The automatic option uses a typical distribution for
                                      proteins and equiprobable distribution for everything else.
                                   2) CG percentage.
                                   3) Explicit distribution of bases, amino acids, or codons
                                      (format: "{'A':10, 'C':40, 'G':40, 'T':10}").
                                   4) Tab-delimited file containing each symbol and its expected
                                      composition in each line.

  -w, --weight NUMBER          The weight of prior data.  Default depends on alphabet length

  -N, --second-data FILENAME   Second dataset file to compare with the input. HeatLogo measures
                               statistical significance at each column by comparing the number of
                               symbols appered in input and second datasets.
                               Accepted formats: clustal, fasta, plain, msf, genbank, nbrf, nexus,
                               phylip, stockholm, intelligenetics, table, array.

  -B, --second-composition COMP
                               The expected composition of the second sample sequences. The supported
                               formats are the same as those of '--composition'.

  -W, --second-weight NUMBER   The prior weight for the second dataset. Default depends on alphabet length.

  -U, --units NUMBER           A unit of entropy ('bits' (default), 'nats', 'digits'), or a unit of free
                               energy ('kT', 'kJ/mol', 'kcal/mol'), or 'probability' for probabilities

Heat Map Options
~~~~~~~~~~~~~~~~

These options affect the format of the sequence logo in the heat map
mode (See also Color Options).

  -S, --stats-test TEST        Specify a statistical test to calculate p-values: 'binomial', 't-test',
                               'z-score' (default: t-test)

  -H, --heat-scheme SCHEME     Specify a heatmap color scheme: 'magenta-cyan', 'red-blue', 'red-green',
                               'thermography' (See also: '--p-color').

      --p-color COLOR PVALUE   Specify a color for a p-value (p-value range: -1 < p < 1).
                               (e.g. --p-color blue 0.5 --p-color red 0.25)

      --hide-colorkey          Toggle switch to hide a heatmap color key on the right side of a
                               sequence logo (default: off)

Color Options
~~~~~~~~~~~~~

Colors can be specified using CSS2 syntax. e.g. 'red', '#FF0000', etc
(See also Heat Map Options).

  -C, --color-scheme SCHEME    Specify a standard color scheme (auto, base pairing, charge, chemistry,
                               classic, dna codon, hydrophobicity, monochrome, rna codon; Heatmap: off)

      --c-color COLOR SYMBOLS DESCRIPTION 
                               Specify symbol colors (e.g. --c-color black AG 'Purine'
                                --c-color red TC 'Pyrimidine'; Heatmap: off)

      --default-color COLOR    Symbol color if not otherwise specified.

Transformations
~~~~~~~~~~~~~~~

Optional transformations of the sequence data.

      --ignore-lower-case      Disregard lower case letters and only count upper case letters in sequences.

      --reverse                Simply reverse input sequences.

      --complement             Complement DNA sequences.

Logo Format Options
~~~~~~~~~~~~~~~~~~~

These options control the format and display of the sequence logo.

  -I, --first-index INDEX      Index of first position in sequence data (default: 1)

  -s, --start INDEX            Lower bound (i.e. start position) of sequence to display

  -e, --end INDEX              Upper bound (i.e. end position) of sequence to display

  -L, --size LOGOSIZE          Specify a standard logo size (small, medium (default), large)

  -n, --stacks-per-line COUNT  Maximum number of logo stacks per logo line (default: 40)

  -t, --title TEXT             Logo title text.

  -l, --label TEXT             A figure label (e.g. '2a').

  -X, --show-xaxis YES/NO      Display sequence numbers along x-axis? (default: True)

  -x, --xlabel TEXT            X-axis label

      --annotate TEXT          A comma separated list of custom stack annotations (e.g. '1,3,4,5,6,7').
                               Annotation list must be same length as sequences.

  -M, --yaxis UNIT             Height of yaxis in units. (Default: Maximum value with uninformative prior.)

  -Y, --show-yaxis YES/NO      Display entropy scale along y-axis? (default: True)

  -y, --ylabel TEXT            Y-axis label (default depends on plot type and units)

  -E, --show-ends YES/NO       Label the ends of the sequence? (default: False)

  -P, --fineprint TEXT         The fine print (default: heatlogo version)

      --ticmarks NUMBER        Distance between ticmarks (default: 1.0)

      --errorbars YES/NO       Display error bars? (default: True)

      --reverse-stacks YES/NO  Draw stacks with largest letters on top? (default: True)

Advanced Format Options
~~~~~~~~~~~~~~~~~~~~~~~

These options provide fine control over the display of the sequence logo.

  -k, --stack-width POINTS     Width of a logo stack (default: 10.8)

      --aspect-ratio POINTS    Ratio of stack height to width (default: 5)

      --box YES/NO             Draw boxes around symbols? (default: no)

      --resolution DPI         Bitmap resolution in dots per inch (DPI) (Default: 96 DPI, except
                               png_print, 600 DPI). Low resolution bitmaps (DPI<300) are antialiased.

      --scale-width YES/NO     Scale the visible stack width by the fraction of symbols in the column?
                               (i.e. columns with many gaps of unknowns are narrow.) (Default: yes)

HeatLogo Server
~~~~~~~~~~~~~~~

Run a standalone webserver on a local port.

      --serve                  Start a standalone HeatLogo server for creating sequence logos.

      --port PORT              Listen to this local port. (Default: 8080)

Examples
--------

For help on the command line interface, type::

    $ ./heatlogo --help

To build a simple logo, execute::

    $ ./heatlogo < seqs.fa > logo.eps

If the input data is codon, type::

    $ ./heatlogo --sequence-type 'codon' < seqs.fa > logo.eps

HeatLogo takes into account expected composition of sequences using the following commands.

- String (for DNA, RNA, protein, and codon)::

    $ ./heatlogo --composition ['auto'| 'equiprobable'|'none'] < seqs.fa > logo.eps

- Explicit distribution of symbols (for DNA, RNA, protein, and codon)::

    $ ./heatlogo --composition "{'A':30, 'C':20, 'G':20, 'T':30}" < seqs.fa > logo.eps

- CG percentage (only applicable for DNA or RNA)::

    $ ./heatlogo --composition 40% < seqs.fa > logo.eps

- Tab-delimited file containing each symbol and its composition (for DNA, RNA, protein, and codon)::

    $ ./heatlogo --composition COMPOSITON_FILENAME < seqs.fa > logo.eps

  * Note: tab-delimited file format

    +--------------+---------------+-------------------+
    |  DNA or RNA  |  protein      |       codon       |
    +==============+===============+===================+
    |   A   24.3   |  A   13.8     |  AAA   0.080      |
    +--------------+---------------+-------------------+
    |   T   24.3   |  C   12.1     |  AAC   0.018      |
    +--------------+---------------+-------------------+
    |   G   25.7   |  D    3.4     |  AAG   0.409      |
    +--------------+---------------+-------------------+
    |   C   25.7   |  E    8.9 ... |  AAT   0.004 ...  |
    +--------------+---------------+-------------------+
    |   (4 lines)  | (20 lines)    |    (64 lines)     |
    +--------------+---------------+-------------------+

For calculating statistical significance of each amino acid,
specify a statistical test as follows::

    $ ./heatlogo --stats-test 'binomial' < seqs.fa > logo.eps

If there is negative or paired dataset against the input, HeatLogo computes P-values
by comparing each column of the input with the column of the negative or paired dataset::

    $ ./heatlogo --second-data FILENAME < seqs.fa > logo.eps

Heatmap schemes implemented in HeatLogo can be changed with the following option::

    $ ./heatlogo --heat-scheme 'red-blue' < seq.fa > logo.eps

To generate position weight matrices of P-values and log-odds ratios, type::

    $ ./heatlogo --pwm-pval FILENAME --pwm-prob FILENAME < seqs.fa > logo.eps

HeatLogo also generates sequence logos with usual color schemes by turning off
the heatmap mode::

    $ ./heatlogo --color-scheme 'base pairing' < seqs.fa > logo.eps

Color Schemes
-------------

Color schemes (i.e. color schemes which are not heatmap mode) of usual
sequence logos are listed as follows.

Amino Acids
~~~~~~~~~~~

All color schemes for amino acids are based on their physicochemical properties.

- Hydrophobicity

+------------------------+------------------------+----------+
| Physicochemical prop.  |       Amino acid       |  Colors  |
+========================+========================+==========+
|      Hydrophilic       |    R, K, D, E, N, Q    |  blue    |
+------------------------+------------------------+----------+
|        Neutral         |    S, G, H, T, A, P    |  green   |
+------------------------+------------------------+----------+
|      Hydrophobic       | Y, V, M, C, L, F, I, W |  black   |
+------------------------+------------------------+----------+

- Chemistry

+------------------------+------------------------+----------+
| Physicochemical prop.  |       Amino acid       |  Colors  |
+========================+========================+==========+
|         Polar          |     G, S, T, Y, C      |  green   |
+------------------------+------------------------+----------+
|        Neutral         |          N, Q          |  purple  |
+------------------------+------------------------+----------+
|         Basic          |        K, R, H         |  blue    |
+------------------------+------------------------+----------+
|        Acidic          |          D, E          |  red     |
+------------------------+------------------------+----------+
|      Hydrophobic       | P, A, W, F, L, I, M, V |  black   |
+------------------------+------------------------+----------+

- Charge

+------------------------+--------------+----------+
| Physicochemical prop.  |  Amino acid  |  Colors  |
+========================+==============+==========+
|        Positive        |   K, R, H    |   blue   | 
+------------------------+--------------+----------+
|        Negative        |     D, E     |   red    |
+------------------------+--------------+----------+

Codons
~~~~~~

The color scheme is based on physicochemical properties of amino acids
encoded by codons according to the standard genetic code table.

+------------------------+--------------+---------------------------+ 
| Physicochemical prop.  |  Amino acid  |  Colors                   | 
+========================+==============+===========================+
|     Polar positive     |    H, K, R   |  light blue               |
+------------------------+--------------+---------------------------+ 
|     Polar negative     |     D, E     |  red                      |
+------------------------+--------------+---------------------------+ 
|     Polar neutral      |  S, T, N, Q  |  green                    |
+------------------------+--------------+---------------------------+
|   Non-polar aliphatic  |  A, V, L, I  |  blue                     |
+------------------------+--------------+---------------------------+
|   Non-polar aromatic   |    F, W, Y   |  magenta                  |
+------------------------+--------------+---------------------------+
|    Special residue     |    P, G, C   |  brown (P, G), yellow (C) | 
+------------------------+--------------+---------------------------+
|       Stop codons      |     ``*``    |  black                    |
+------------------------+--------------+---------------------------+

