#!/usr/bin/env python

# -------------------------------- WebLogo --------------------------------
#
#  Copyright (c) 2003-2004 The Regents of the University of California.
#  Copyright (c) 2005 Gavin E. Crooks
#  Copyright (c) 2006-2011, The Regents of the University of California, through 
#  Lawrence Berkeley National Laboratory (subject to receipt of any required
#  approvals from the U.S. Dept. of Energy).  All rights reserved.
#
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

# WebLogo Command Line Interface

import sys

from color import *
from colorscheme import ColorScheme, ColorGroup
from colorscheme import HeatScheme, HeatGroup
from corebio.utils import *
from string import Template

from corebio import seq_io
from corebio.seq import Seq, SeqList

import os 
from corebio.utils.deoptparse import DeOptionParser
from optparse import OptionGroup

from weblogolib import LogoOptions, LogoData, LogoFormat, writePssm
from weblogolib import parse_prior, description, release_description, formatters, default_formatter
from weblogolib import codon_dna_alphabet, codon_rna_alphabet
from weblogolib import std_alphabets, std_units, std_sizes, std_color_schemes
from weblogolib import stats_tests, heat_color_schemes
from weblogolib import read_seq_data

# ====================== Main: Parse Command line =============================
def main(): 
    """HeatLogo command line interface """
    
    # ------ Parse Command line ------
    parser = _build_option_parser()                
    (opts, args) = parser.parse_args(sys.argv[1:])
    if args : parser.error("Unparsable arguments: %s " % args)
    
    if opts.serve:
        httpd_serve_forever(opts.port) # Never returns?
        sys.exit(0) 
            
    # ------ Create Logo ------
    try:
        data = _build_logodata(opts)
        format = _build_logoformat(data, opts)
        
        formatter = opts.formatter
        formatter(data, format, opts.fout)

        if opts.pwm_prob:
            writePssm(data, format, opts.pwm_prob)
        if opts.pwm_pval:
            writePssm(data, format, opts.pwm_pval, True)

    except ValueError, err :
        print >>sys.stderr, 'Error:', err
        sys.exit(2)
    except KeyboardInterrupt, err:
        sys.exit(0)
# End main()            
        

def httpd_serve_forever(port=8080) :
    """ Start a webserver on a local port."""
    import BaseHTTPServer
    import CGIHTTPServer 
    
    class __HTTPRequestHandler(CGIHTTPServer.CGIHTTPRequestHandler):
        # Modify CGIHTTPRequestHandler so that it will run the cgi script directly, instead of exec'ing
        # This bypasses the need for the cgi script to have execute permissions set,
        # since distutils install does not preserve permissions.
        def is_cgi(self) :
            self.have_fork = False          # Prevent CGIHTTPRequestHandler from using fork
            if self.path == "/create.cgi": 
                self.cgi_info = '', 'create.cgi'
                return True
            return False
        def is_python(self,path):           # Let CGIHTTPRequestHandler know that cgi script is python
            return True
            
    
    # Add current directory to PYTHONPATH. This is
    # so that we can run the standalone server
    # without having to run the install script.      
    pythonpath = os.getenv("PYTHONPATH", '')
    pythonpath += ":" + os.path.abspath(sys.path[0]).split()[0]
    os.environ["PYTHONPATH"] = pythonpath

    htdocs = resource_filename(__name__, 'htdocs', __file__)
    os.chdir(htdocs) 

    HandlerClass = __HTTPRequestHandler
    ServerClass = BaseHTTPServer.HTTPServer
    httpd = ServerClass(('', port), HandlerClass)
    print "HeatLogo server running at http://localhost:%d/" % port
    
    try :
        httpd.serve_forever()
    except KeyboardInterrupt:
        sys.exit(0)
# end httpd_serve_forever()  


def _build_logodata_core(fin, fin_compos, fin_weight, options, second_data=None):
    motif_flag=False
    isCodon = False
   
    try:
        # Try reading data in transfac format first.
        from corebio.matrix import Motif
        motif = Motif.read_transfac(fin, alphabet=options.alphabet)
        motif_flag = True
    except ValueError, motif_err:
        # Failed reading Motif, try reading as multiple sequence data.
        if options.alphabet is True:
            isCodon = True
            options.alphabet = None
        seqs = read_seq_data(fin,
            options.input_parser.read,
            alphabet=options.alphabet,
            ignore_lower_case = options.ignore_lower_case)

    if motif_flag :
        if options.ignore_lower_case:
            raise ValueError("option --ignore-lower-case incompatible with matrix input")
        if options.reverse: motif.reverse()
        if options.complement: motif.complement()

        if not isCodon:
            prior, compos = parse_prior(fin_compos, motif.alphabet, fin_weight)
            data = LogoData.from_counts(motif.alphabet, motif, options.stats_func, prior, compos, second_data)
        else:
            raise ValueError("option --sequence-type 'codon' incompatible with matrix input")
    else :
        if options.codon_frame < 0 and isCodon:
            options.reverse = True
            options.complement = True

        if options.reverse: 
            seqs = SeqList( [s.reverse() for s in seqs], seqs.alphabet )
        
        if options.complement:
            seqs = SeqList( [Seq(s,seqs.alphabet).complement() for s in seqs], seqs.alphabet )

        if isCodon:
            if abs(options.codon_frame) > 1:
                beg = abs(options.codon_frame)-1
                end = beg + int( (len(seqs[0])-beg)/3 ) * 3
                seqs = SeqList([s[beg:end] for s in seqs], seqs.alphabet)
            if std_alphabets['dna'] == seqs.alphabet:
                seqs.alphabet = codon_dna_alphabet
            elif std_alphabets['rna'] == seqs.alphabet:
                seqs.alphabet = codon_rna_alphabet

        prior,compos = parse_prior(fin_compos, seqs.alphabet, fin_weight)
        data = LogoData.from_seqs(seqs, options.stats_func, prior, compos, second_data)
    return data

  
def _build_logodata(options):
    if options.second_data is not None:
        ngdata = _build_logodata_core(options.second_data,
                                      options.second_composition,
                                      options.second_weight, options)
    fin = options.fin
    if fin is None:
        from StringIO import StringIO
        fin = StringIO(sys.stdin.read())
    if options.second_data is not None:
        data = _build_logodata_core(fin, options.composition, options.weight, options, ngdata)
    else:
        data = _build_logodata_core(fin, options.composition, options.weight, options)
    return data

             
def _build_logoformat(logodata, opts) :
    """ Extract and process relevant option values and return a 
    LogoFormat object.""" 

    args = {}  
    direct_from_opts = [
        "stacks_per_line", 
        "logo_title",
        "yaxis_label", 
        "show_xaxis",
        "show_yaxis",        
        "xaxis_label", 
        "show_ends",
        "fineprint",  
        "show_errorbars", 
        "show_boxes",  
        "yaxis_tic_interval", 
        "resolution",      
        "alphabet",
        "show_ends",
        "default_color",
        "color_scheme",
        "heat_scheme",
        "p_colors",
        "stats_func",
        "show_colorkey",
        "unit_name",
        "logo_label",
        "yaxis_scale",
        "first_index",
        "logo_start",
        "logo_end",
        "scale_width", 
        "annotate",
        "stack_width",
        "stack_aspect_ratio",
        "reverse_stacks"
        ]
  
    for k in direct_from_opts:
        args[k] = opts.__dict__[k]

#    logo_size = copy.copy(opts.__dict__['logo_size'])
#    size_from_opts = ["stack_width", "stack_height"]
#    for k in size_from_opts :
#        length = getattr(opts, k)
#        if length : setattr( logo_size, k, length )
#   args["size"] = logo_size    

    if opts.colors:
        color_scheme = ColorScheme()
        for color, symbols, desc in opts.colors:
            try :
                #c = Color.from_string(color)
                color_scheme.groups.append( ColorGroup(symbols, color, desc)  )
            except ValueError : 
                raise ValueError("error: option --ch-color: invalid value: '%s'" % color )
        args["color_scheme"] = color_scheme

    if opts.p_colors:
        heat_scheme = HeatScheme()
        for color, pvalues in opts.p_colors:
            try :
                heat_scheme.groups.append( HeatGroup(pvalues, color) )
            except ValueError :
                raise ValueError("error: option --p-color: invalid value: '%s'" % color )
        args["heat_scheme"] = heat_scheme
  
    if opts.annotate:
        args["annotate"] = opts.annotate.split(',')
    
    logooptions = LogoOptions()
    for a, v in args.iteritems() :
        setattr(logooptions,a,v)
    
    theformat =  LogoFormat(logodata, logooptions)
    return theformat

    
# ========================== OPTIONS ==========================
def _build_option_parser() :
    defaults = LogoOptions()
    parser = DeOptionParser(usage="%prog [options] < sequence_data.fa > sequence_logo.eps",
        description = description,
        version     = release_description,
        add_verbose_options = False
        )
        
    io_grp = OptionGroup(parser, "Input/Output Options",)
    data_grp = OptionGroup(parser, "Logo Data Options",)
    heat_grp = OptionGroup(parser, "Heat Map Options",
        "These options affect the format of the sequence logo in the heat map mode (See also Color Options).")
    color_grp = OptionGroup(parser, "Color Options",
        "Colors can be specified using CSS2 syntax. e.g. 'red', '#FF0000', etc (See also Heat Map Options).")
    trans_grp = OptionGroup(parser, "Transformations",
        "Optional transformations of the sequence data.")
    format_grp = OptionGroup(parser, "Logo Format Options",
        "These options control the format and display of the sequence logo.")
    advanced_grp = OptionGroup(parser, "Advanced Format Options", 
        "These options provide fine control over the display of the sequence logo.")
    server_grp = OptionGroup(parser, "HeatLogo Server",
        "Run a standalone webserver on a local port.")

    parser.add_option_group(io_grp)
    parser.add_option_group(data_grp)
    parser.add_option_group(heat_grp)
    parser.add_option_group(color_grp)
    parser.add_option_group(trans_grp)      
    parser.add_option_group(format_grp)  
    parser.add_option_group(advanced_grp)
    parser.add_option_group(server_grp)
    
    # ========================== IO OPTIONS ==========================
    
    io_grp.add_option("-i", "--input",
        dest="fin",
        action="store",
        type="file_in",
        default=None,
        help="Sequence input file (default: stdin)",
        metavar="FILENAME")

    # Add position weight matrix formats to input parsers by hand
    fin_choices = dict(seq_io.format_names())
    fin_choices['transfac'] = 'transfac'
    formatters_string = formatters.keys()
    io_grp.add_option("-f", "--input-format", 
        dest="input_parser",
        action="store", type ="dict",
        default = seq_io,
        choices = fin_choices,       # seq_io.format_names(),
        help="Type of multiple sequence alignment or position weight matrix file: (%s, transfac)" % 
           ', '.join([ f.names[0] for f in seq_io.formats]),
        metavar="FORMAT")

    io_grp.add_option("-o", "--output",
        dest="fout",
        type="file_out",
        default=sys.stdout,
        help="Output file (default: stdout)",
        metavar="FILENAME")

    io_grp.add_option("-F", "--output-format",
        dest="formatter",
        action="store",
        type="dict",
        metavar= "FORMAT",
        choices=formatters,
        help="Format of output: eps (default), png, png_print, pdf, jpeg, svg, logodata",
        default = default_formatter)

    io_grp.add_option("", "--pwm-prob",
        dest="pwm_prob",
        action="store",
        type="string",
        metavar="FILENAME",
        default=None,
        help="Output position weight matrix composed of logarithmic probabilities")

    io_grp.add_option("", "--pwm-pval",
        dest="pwm_pval",
        action="store",
        type="string",
        metavar="FILENAME",
        default=None,
        help="Output position weight matrix of P-values")


    # ========================== Data OPTIONS ==========================
 
    data_grp.add_option("-A", "--sequence-type",
        dest="alphabet",
        action="store",
        type="dict",
        choices = std_alphabets,
        help="The type of sequence data: 'protein', 'rna', 'dna', or 'codon' (See also: --codon-frame).",
        metavar="TYPE")

                        
    data_grp.add_option("-a", "--alphabet",
        dest="alphabet",
        action="store",
        help="The set of symbols to count, e.g. 'AGTC'. All characters not in the alphabet are ignored. If neither the alphabet nor sequence-type are specified then heatlogo will examine the input data and make an educated guess. (See also: --sequence-type, --ignore-lower-case)" )

    data_grp.add_option("-R", "--codon-frame",
        dest="codon_frame",
        action="store",
        type="int",
        default=defaults.codon_frame,
        help="Codon reading frame (default: +1): [+1, +2, +3] indicate the forward strand frame, and [-1, -2, -3] indicate the reverse strand frame",
        metavar="NUMBER")

    data_grp.add_option("-b", "--composition",
        dest="composition",
        action="store",
        type="string",
        default = "auto",
        help="The expected composition of the sequences. HeatLogo accepts four types of formats: (1) String: 'auto' (default), 'equiprobable', or 'none' (do not perform any compositional adjustment, Heatmap: off). The automatic option uses a typical distribution for proteins and equiprobable distribution for everything else. (2) CG percentage. (3) Explicit distribution of bases, amino acids, or codons (format: \"{'A':10, 'C':40, 'G':40, 'T':10}\"). (4) Tab-delimited file containing each symbol and its expected composition in each line.",
        metavar="COMP.")

    data_grp.add_option("-w", "--weight",
        dest="weight",
        action="store",
        type="float",
        default = None,
        help="The weight of prior data.  Default depends on alphabet length",
        metavar="NUMBER")

    data_grp.add_option("-N", "--second-data",
        dest="second_data",
        action="store",
        type="file_in",
        default=None,
        help="Second dataset file to compare with the input. HeatLogo measures statistical significance at each column by comparing the number of symbols appered in input and second datasets. Accepted formats: %s." % ', '.join([ f.names[0] for f in seq_io.formats]),
        metavar="FILENAME")

    data_grp.add_option("-B", "--second-composition",
        dest="second_composition",
        action="store",
        type="string",
        default="auto",
        help="The expected composition of the second sample sequences. The supported formats are the same as those of '--composition'.",
        metavar="2nd-COMP.")

    data_grp.add_option("-W", "--second-weight",
        dest="second_weight",
        action="store",
        type="float",
        default = None,
        help="The prior weight for the second dataset.  Default depends on alphabet length",
        metavar="NUMBER")

    data_grp.add_option("-U", "--units",
        dest="unit_name",
        action="store",
        choices = std_units.keys(),
        type="choice",
        default = defaults.unit_name,
        help="A unit of entropy ('bits' (default), 'nats', 'digits'), or a unit of free energy ('kT', 'kJ/mol', 'kcal/mol'), or 'probability' for probabilities",
        metavar = "NUMBER")
       
    # ======================== Heat Map OPTIONS =========================

    stats_choices = stats_tests.keys()
    stats_choices.sort()
    heat_grp.add_option("-S", "--stats-test",
        dest="stats_func",
        action="store",
        type="dict",
        choices = stats_tests,
        metavar="TEST",
        default = defaults.stats_func,
        help="Specify a statistical test to calculate p-values: %s (default: t-test)" % \
             ", ".join(stats_choices) )

    heat_color_choices = heat_color_schemes.keys()
    heat_color_choices.sort()
    heat_grp.add_option("-H", "--heat-scheme",
        dest="heat_scheme",
        action="store",
        type="dict",
        choices = heat_color_schemes,
        metavar="SCHEME",
        default = defaults.heat_scheme,
        help="Specify a heatmap color scheme (%s) (See also: '--p-color')" % \
             ", ".join(heat_color_choices) )

    heat_grp.add_option("", "--p-color",
        dest="p_colors",
        action="append",
        metavar="COLOR PVALUE",
        nargs = 2,
        default=[],
        help="Specify a color for a p-value (p-value range: -1 < p < 1), e.g. --p-color blue 0.5 --p-color red 0.25")

    heat_grp.add_option("", "--hide-colorkey",
        dest="show_colorkey",
        action="store_false",
        default=True,
        help="Toggle switch to hide a heatmap color key on the right side of a sequence logo (default: off)" )

    # ========================== Color OPTIONS ==========================

    color_scheme_choices = std_color_schemes.keys()
    color_scheme_choices.sort()    
    color_grp.add_option("-C", "--color-scheme",
        dest="color_scheme",
        action="store",
        type ="dict",
        choices = std_color_schemes,
        metavar = "SCHEME",
        default = None, # Auto
        help="Specify a standard color scheme (%s) (Heatmap: off)" % \
             ", ".join(color_scheme_choices) )
            
    color_grp.add_option("", "--c-color",
        dest="colors",
        action="append",
        metavar="COLOR SYMBOLS DESCRIPTION ",
        nargs = 3,
        default=[],
        help="Specify symbol colors, e.g.  --c-color black AG 'Purine' --c-color red TC 'Pyrimidine', (Heatmap: off)")

    color_grp.add_option("", "--default-color",
        dest="default_color",
        action="store",
        metavar="COLOR",
        default= defaults.default_color,
        help="Symbol color if not otherwise specified.")

    # ========================== Transformation OPTIONS ==========================
    
    # FIXME Add test?
    trans_grp.add_option("", "--ignore-lower-case",
        dest="ignore_lower_case",
        action="store_true",
        default=False,
        help="Disregard lower case letters and only count upper case letters in sequences."
       )

    trans_grp.add_option("", "--reverse",
        dest="reverse",
        action="store_true",
        default=False,
        help="reverse sequences",
        )

    trans_grp.add_option("", "--complement",
        dest="complement",
        action="store_true",
        default=False,
        help="complement DNA sequences" )


    # ========================== FORMAT OPTIONS ==========================

    format_grp.add_option("-I", "--first-index",
        dest="first_index",
        action="store",
        type="int",
        default = 1,
        help="Index of first position in sequence data (default: 1)",
        metavar="INDEX")

    format_grp.add_option("-s", "--start",
        dest="logo_start",
        action="store",
        type="int",
        help="Lower bound (i.e. start position) of sequence to display",
        metavar="INDEX")
    
    format_grp.add_option("-e", "--end",
        dest="logo_end",
        action="store",
        type="int",
        help="Upper bound (i.e. end position) of sequence to display",
        metavar="INDEX")

    format_grp.add_option("-L", "--size",
        dest="stack_width",
        action="store",
        type ="dict",
        choices = std_sizes,
        metavar = "LOGOSIZE",
        default = defaults.stack_width,
        help="Specify a standard logo size (small, medium (default), large)" )

    format_grp.add_option("-n", "--stacks-per-line",
        dest="stacks_per_line",
        action="store",
        type="int",
        help="Maximum number of logo stacks per logo line. (default: %default)",
        default = defaults.stacks_per_line,
        metavar="COUNT")

    format_grp.add_option("-t", "--title",
        dest="logo_title",
        action="store",
        type="string",
        help="Logo title text.",
        default = defaults.logo_title,
        metavar="TEXT")

    format_grp.add_option("-l", "--label",
        dest="logo_label",
        action="store",
        type="string",
        help="A figure label, e.g. '2a'",
        default = defaults.logo_label,
        metavar="TEXT")

    format_grp.add_option("-X", "--show-xaxis",
        action="store",
        type = "boolean",
        default= defaults.show_xaxis,
        metavar = "YES/NO",
        help="Display sequence numbers along x-axis? (default: %default)")
                       
    format_grp.add_option("-x", "--xlabel",
        dest="xaxis_label",
        action="store",
        type="string",
        default = defaults.xaxis_label,
        help="X-axis label",
        metavar="TEXT")

    format_grp.add_option("", "--annotate",
            dest="annotate",
            action="store",
            type="string",
            default = None,
            help="A comma separated list of custom stack annotations, e.g. '1,3,4,5,6,7'.  Annotation list must be same length as sequences.",
            metavar="TEXT")

    format_grp.add_option("-M", "--yaxis",
        dest="yaxis_scale",
        action="store",
        type="float",
        help="Height of yaxis in units. (Default: Maximum value with uninformative prior.)",
        metavar = "UNIT") 

    format_grp.add_option("-Y", "--show-yaxis",
        action="store",
        type = "boolean",
        dest = "show_yaxis",
        default= defaults.show_yaxis,
        metavar = "YES/NO",
        help="Display entropy scale along y-axis? (default: %default)")

    format_grp.add_option("-y", "--ylabel",
        dest="yaxis_label",
        action="store",
        type="string",
        help="Y-axis label (default depends on plot type and units)",
        metavar="TEXT")

    format_grp.add_option("-E", "--show-ends",
        action="store",
        type = "boolean",
        default= defaults.show_ends,
        metavar = "YES/NO",
        help="Label the ends of the sequence? (default: %default)")
        
    format_grp.add_option("-P", "--fineprint",
        dest="fineprint",
        action="store",
        type="string",
        default= defaults.fineprint,
        help="The fine print (default: heatlogo version)",
        metavar="TEXT")

    format_grp.add_option("", "--ticmarks",
        dest="yaxis_tic_interval",
        action="store",
        type="float",
        default= defaults.yaxis_tic_interval,
        help="Distance between ticmarks (default: %default)",
        metavar = "NUMBER")

        
    format_grp.add_option("", "--errorbars",
        dest = "show_errorbars",
        action="store",
        type = "boolean",
        default= defaults.show_errorbars,
        metavar = "YES/NO",
        help="Display error bars? (default: %default)")
     
    format_grp.add_option("", "--reverse-stacks",
        dest = "reverse_stacks",
        action="store",
        type = "boolean",
        default= defaults.show_errorbars,
        metavar = "YES/NO",
        help="Draw stacks with largest letters on top? (default: %default)")
 

    # ========================== Advanced options =========================   
                
    advanced_grp.add_option("-k", "--stack-width",
        dest="stack_width",
        action="store",
        type="float",
        default= defaults.stack_width,
        help="Width of a logo stack (default: %s)"% defaults.stack_width,
        metavar="POINTS" )

    advanced_grp.add_option("", "--aspect-ratio",
        dest="stack_aspect_ratio",
        action="store",
        type="float",
        default= defaults.stack_aspect_ratio ,
        help="Ratio of stack height to width (default: %s)"%defaults.stack_aspect_ratio,
        metavar="POINTS" )    

    advanced_grp.add_option("", "--box",
        dest="show_boxes",
        action="store",
        type = "boolean",
        default=False, 
        metavar = "YES/NO",
        help="Draw boxes around symbols? (default: no)")

    advanced_grp.add_option("", "--resolution",
        dest="resolution",
        action="store",
        type="float",
        default=96,
        help="Bitmap resolution in dots per inch (DPI).  (Default: 96 DPI, except png_print, 600 DPI) Low resolution bitmaps (DPI<300) are antialiased.",
        metavar="DPI")  

    advanced_grp.add_option("", "--scale-width",
        dest="scale_width",
        action="store",
        type = "boolean",
        default= True, 
        metavar = "YES/NO",
        help="Scale the visible stack width by the fraction of symbols in the column?  (i.e. columns with many gaps of unknowns are narrow.)  (Default: yes)")
   
    # ========================== Server options =========================   
    server_grp.add_option("", "--serve",
        dest="serve",
        action="store_true",
        default= False,
        help="Start a standalone HeatLogo server for creating sequence logos.")    
    
    server_grp.add_option("", "--port",
        dest="port",
        action="store",
        type="int",
        default= 8080,
        help="Listen to this local port. (Default: %default)",
        metavar="PORT")

    return parser
    
    # END _build_option_parser

##############################################################
