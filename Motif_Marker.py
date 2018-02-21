#!/Users/Adrian/miniconda3/bin/python3

## Motif Marker:
## Adrian Bubie
## 2/11/18
## --------------------
## This program detects motifs (selected sequences) of interest from sequences provided in a fasta format. 
## As the intended use case is to identify motifs around exon boundries that are potential regulators for splicing events, 
## this program maps detected motifs onto a representation of the intron/exons sequences provided in the input file, and
## exports these maps as an SVG graphic.

import textwrap
import random
import re
import argparse as ap
import cairo
import math

def get_arguments():
    parser = ap.ArgumentParser(prog='./Motif_Marker.py', formatter_class=ap.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
    Motif Sequence Marker
    ---------------------
    Locates and counts specified motif sequences from given gene sequences
    in FASTA format, then maps motifs' locations in relation to Intron/Exon of given sequence.
    
    Takes input of FASTA file with any exonic sequences capitalized, and intronic sequences lowercase.
    Requires motifs to search for and size of intronic sequence flanking exons to be searched to be specified by the user. 
    '''))
    parser.add_argument("-f", help="Fasta file to be processed. Exon sequences must be capitalized, and gene names contained in sequence ID. Must include absolute path the file. <str>", required=True, type=str)
    parser.add_argument("-w", help="Intronic flanking sequence window sized (applies to both sides of exonic regions). Default is 200bp. <int>", required=False, type=int, default=200)
    parser.add_argument("-m", help="Motif file to define motif sequences to query for (one Motif per line). Motifs must use *IUPAC* Nucleotide codes. Must include absolute path to file <str>", required=True, type=str)
    parser.add_argument("-s", help="Optional summary output file containing the count of each motif per fasta sequence. File is created in location of imput file if set to 'True' <str>", required=False, type=bool, default=False)
    #parser.add_argument("-colors", help="Set color palate to use for motif markers (note: whole sequences will always be represented in black). Provide as list of hex codes (#)", required=False, type=str, default='')
    return parser.parse_args()


class fasta_sequence():
    '''Fasta Sequence object: collection of sequence lines for a single fasta ID.
       Object collects lines into single string for motif searching, collects seq length, and more.'''
    
    def __init__(self, lines):
        self.gene = lines[0].split(' ')[0].strip('>')
        self.seq = ''.join([line.strip('\n') for line in lines if line.startswith('>') != True])
    
    def tot_seq_len(self):
        return int(len(self.seq))
    
    def exon_bounds(self):
        exon_st = self.seq.find(re.search('[ATCG]+',self.seq)[0])
        exon_ed = (self.seq[exon_st:].find(re.search('[atcg]+',self.seq[exon_st:])[0]))+ exon_st
        return [exon_st, exon_ed]

 
def iupac_interp(motifs_file):
    '''IUPAC Interpreter: Takes motif file, parses out motifs, and returns list of regex search terms for each motif
       based on IUPAC nomencalture.'''
    
    # Keys of dict are nucleotide code, and values are bases:
    iupac = {"A":'[Aa]',"C":'[Cc]',"G":'[Gg]',"T":'[TUtu]',"U":'[TUtu]',"R":'[AGag]',"Y":'[CTct]',"S":'[GCgc]',"W":'[ATat]',"K":'[GTgt]',"M":'[ACac]',"B":'[CGTcgt]',"D":'[AGTagt]',"H":'[ACTact]',"V":'[ACGacg]',"N":'[A-Za-z]',
             "a":'[Aa]',"c":'[Cc]',"g":'[Gg]',"t":'[TUtu]',"u":'[TUtu]',"r":'[AGag]',"y":'[CTct]',"s":'[GCgc]',"w":'[ATat]',"k":'[GTgt]',"m":'[ACac]',"b":'[CGTcgt]',"d":'[AGTagt]',"h":'[ACTact]',"v":'[ACGacg]',"n":'[A-Za-z]'}
    
    with open(motifs_file, 'r') as mtfs:
        search_terms = []       # Create a list to store the returned search terms
        line = mtfs.readline()  
        while line:             # For each motif in the file
            st = ''
            mt = str(line).strip('\n')
            for char in mt:             # For each character in the motif
                if char in iupac.keys():
                    st = st+iupac[char]   # If the character is in the IUPAC dict, add it to the current search term
                else:
                    raise ValueError('Error: motif contains character not in IUPAC nucleotide codes; motif cannot be translated') # Throw exception if the motif contains character not in IUPAC
            
            search_terms.append(st)
            line = mtfs.readline()
    
    return search_terms       # Return the search terms


def window_trim(fasta, window_size):
    '''Window Trimmer: Takes in a fasta gene sequence and trims the flagging intronic sequences surrounding the exonic,
       uppercase sequence to the appropriate window size. If the window size is bigger than the intronic segments, the 
       entire sequence is returned.'''
    
    exon_bounds = fasta.exon_bounds()
    
    if len(fasta.seq[:exon_bounds[0]]) > window_size:
        search_start = len(fasta.seq[:exon_bounds[0]])-window_size
    else:
        search_start = 0
    
    if len(fasta.seq[exon_bounds[1]:]) > window_size:
        search_end = exon_bounds[1]+window_size
    else:
        search_end = fasta.tot_seq_len()
    
    trimmed_seq = fasta.seq[search_start:search_end]
    
    return trimmed_seq


def Motif_search(fasta, window, motifs):
    '''Motif Search: Takes in a fasta gene sequence and list of motifs to be searched for, then searches the trimmed
       version of the sequence for each motif in the set of motifs, and saves the motif positions to a dictionary. Returns
       this dictionary of search results.'''
    
    # Start by getting the window trimmed seq:
    w_seq = window_trim(fasta, window)
    
    # For each motif, search the sequence for the positions of these motifs, and store the positions:
    motif_positions = {}
    for mot in motifs:
        motif_positions[mot]=[x.start(0) for x in re.finditer(mot, w_seq)]

    return motif_positions


def summary_out(results):
    '''Summary File: Takes in the motif search results and produces a text file of the number and location of each motif
       type for each of the searched sequences in the fasta. Summary file is written out to the current directory 
       (where script is executed from)'''
    
    with open('./Motif_Search_Summary.txt','w') as out:
        out.write('## Summary of Motif Search:\n')
        for res in results:
            out.write('Gene: '+res[0].gene+'\n')
            for key in res[1].keys():
                out.write('\tMotif: '+key+'\n')
                out.write('\tNumber found: '+str(len(res[1][key]))+'\n')
                out.write('\tPosition(s) in sequence: '+str(res[1][key])+'\n')
    

def py_draw(search_results, Motifs):
    '''Pycairo Graphing: Takes in the motif search results and draws to-scale visual pertaining to the sequence, exons
       and motifs found. Produces an .svg of the graph in the current directory (where script is executed from)'''
    
    no_of_graphs = len(search_results) # Number of graphs to be drawn
    
    max_size = 0                       # Get the max size among the sequences you have to determine surface width
    for res in search_results:
        if len(res[0].seq) > max_size:
            max_size = len(res[0].seq)
    
    # Get Motif color schemes for legend:
    col = 1
    motif_color = {}
    for motif in Motifs:
        r = 0.2+(col/5)
        g = 0.8-(col/5)
        b = 0.2+(col/10)
        motif_color[motif]=[r,g,b]
        col = col+1 
    
    
    # Start Graph drawing
    surface = cairo.SVGSurface("./exon_graphs.svg", max_size+45, (no_of_graphs*80)+100) # width, height for dimensions
    context = cairo.Context(surface)
    
    for i in range(0,len(search_results)):    # For each sequence with results:
        context.set_line_width(1)
        
        # Add Gene name before graph
        context.set_font_size(10)
        context.set_source_rgb(0, 0, 0)
        context.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
        context.move_to(5,30*(i+1)+50)  # x,y coords; spacing graphs by 30px
        context.show_text(search_results[i][0].gene)
        
        # Draw Line to represent Sequence
        context.move_to(35,30*(i+1)+50) # x,y coords; spacing graphs by 30px
        context.line_to(len(search_results[i][0].seq)+35,30*(i+1)+50) # Drawing line px length of sequence
        context.stroke()
        
        # Draw Rectangle to represent Exon
        exon_coords = search_results[i][0].exon_bounds()  # Get coordinates for the exon in the sequence
        context.rectangle(exon_coords[0]+35,(30*(i+1)+50)-10,exon_coords[1]-exon_coords[0],18)       # x,y top corner, followed by width, height of rec.
        context.fill()
        
        # Add Motifs as small rectangles in appropriate locs
        col = 1
        for key in search_results[i][1].keys():
            context.set_source_rgb(motif_color[key][0],motif_color[key][1],motif_color[key][2])
            mot_width = key.count('[')
            for loc in search_results[i][1][key]:
                context.rectangle(loc+35,(30*(i+1)+50)-10,mot_width,20)
                context.fill()
            col = col+1
    
    # Draw legend
    context.move_to(10,(no_of_graphs*80))
    context.set_source_rgb(0, 0, 0)
    context.show_text("Key")
    context.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    for i in range(0,len(Motifs)):
        context.set_source_rgb(motif_color[Motifs[i]][0],motif_color[Motifs[i]][1],motif_color[Motifs[i]][2])
        context.rectangle(10,(no_of_graphs*80)+((i+1)*15),5,10)
        context.fill()
        
        context.move_to(20,(no_of_graphs*80)+((i+1)*15)+7)
        context.set_source_rgb(0, 0, 0)
        context.show_text(Motifs[i])
        context.move_to(10,(no_of_graphs*80)+((i+1)*15))
        
    surface.finish()


####################
## Main function: ##
####################

# Grab Argument Values:
args = get_arguments()

Fasta = open(args.f,'r')
window = args.w
Motifs = iupac_interp(args.m)

#Motifs = iupac_interp('/Users/Adrian/BGMP/Motif_marker/test/motifs.txt')
#Fasta = open('/Users/Adrian/BGMP/Motif_marker/test/fasta_t1.fa','r')

line = Fasta.readline()
seq_results = []

while line:
    to_process = []
    to_process.append(line)
    line = Fasta.readline()
    while (line.startswith('>') != True) and line != '':
        to_process.append(line)
        line = Fasta.readline()
    
    results = Motif_search(fasta_sequence(to_process), 200, Motifs)
    
    seq_results.append([fasta_sequence(to_process),results])

Fasta.close()


# Create Summary File if flag is True:

if args.s == True:
	summary_out(seq_results)
    
# Call function to draw motif graphs:

py_draw(seq_results, Motifs)

