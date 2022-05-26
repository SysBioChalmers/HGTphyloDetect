#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::TreeIO;

my $infile = $ARGV[0]; #newick tree must end in .tree

$infile =~ /^(.*)\.tree$/;
my $query = $1;
my %color_hash;
popcolhash();

my $outfile1 = "${query}-color.txt";
my $outfile2 = "${query}-font.txt";

open (OFIL1, '>', $outfile1) or die "Couldn't write to file $outfile1: $!\n";
open (OFIL2, '>', $outfile2) or die "Couldn't write to file $outfile2: $!\n";
popcolorfile();
popfontfile();

my $in = new Bio::TreeIO(-file => "$infile", -format => 'newick');
while( my $tree = $in->next_tree ) {
	for my $node ( $tree->get_nodes ) {
		if ( $node->is_Leaf ) {
			my $name = $node->id();

			my $lineage = '';
			if ( $name =~ /.+_(other_[A-Za-z0-9]+)$/ or $name =~ /.+_([A-Za-z0-9]+)$/){
				$lineage = $1;
			}
			
			if ($name =~ /^QUERY/i){
				print OFIL2 "${name},label,#e31a1c,bold,1\n";
				print OFIL1 "${name} \#e31a1c Saccharomycotina\n";
			}else{
				print OFIL2 "${name},label,#000000,normal,1\n";		
			}
			
			if (exists $color_hash{$lineage}) {
			my $color = $color_hash{$lineage};
			print OFIL1 "${name} ${color}\n";
			}else{unless ($name =~ /^QUERY/i) {print OFIL1 "${name} \#969696 Other\n";}}
			
		}
	}
}

close OFIL1;
close OFIL2;

sub popcolorfile
{
print OFIL1 "DATASET_COLORSTRIP
#In colored strip datasets, each ID is associated to a color box/strip and can have an optional label. Color can be specified in hexadecimal, RGB or RGBA notation. When using RGB or RGBA notation, you cannot use COMMA as the dataset separator

#lines starting with a hash are comments and ignored during parsing

#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file (except in the SEPARATOR line, which uses space).

#SEPARATOR TAB
#SEPARATOR COMMA
SEPARATOR SPACE

#label is used in the legend table (can be changed later)
DATASET_LABEL Taxonomy

#dataset color (can be changed later)
COLOR #000000

#=================================================================#
#                    OPTIONAL SETTINGS                            #
#=================================================================#

#If COLOR_BRANCHES is set to 1, branches of the tree will be colored according to the colors of the strips above the leaves.
#When all children of a node have the same color, it will be colored the same, ie. the color will propagate inwards towards the root.
COLOR_BRANCHES 1

#each dataset can have a legend, which is defined below
#for each row in the legend, there should be one shape, color and label
#shape should be a number between 1 and 5:
#1: square
#2: circle
#3: star
#4: right pointing triangle
#5: left pointing triangle

LEGEND_TITLE Taxonomy
LEGEND_SHAPES 1 1 1 1
LEGEND_COLORS #e31a1c #28b1aa #2D33E7 #969696
LEGEND_LABELS Saccharomycotina Fungi Bacteria Other 
# LEGEND_SHAPES 1 1 1 1 1 1 1 1 1 1 1 1
# LEGEND_COLORS #e31a1c #fb9a99 #ff7f00 #28b1aa #b15928 #6a3d9a #cab2d6 #1f78b4 #a6cee3 #33a02c #ffff99 #969696 
# LEGEND_LABELS Leotiomycetes Sordariomycetes Eurotiomycetes Dothideomycetes other_Pezizomycotina other_Ascomycota other_Fungi other_Opisthokonta other_Eukaryota Bacteria Archaea Viruses 

#=================================================================#
#     all other optional settings can be set or changed later     #
#           in the web interface (under 'Datasets' tab)           #
#=================================================================#

#width of the colored strip
STRIP_WIDTH 35

#left margin, used to increase/decrease the spacing to the next dataset. Can be negative, causing datasets to overlap.
#MARGIN 5

#border width; if set above 0, a border of specified width (in pixels) will be drawn around the color strip 
#BORDER_WIDTH 2

#border color; used when BORDER_WIDTH is above 0
#BORDER_COLOR #969696

#always show internal values; if set, values associated to internal nodes will be displayed even if these nodes are not collapsed. It could cause overlapping in the dataset display.
#SHOW_INTERNAL 1

#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages

#=================================================================#
#       Actual data follows after the DATA keyword                #
#=================================================================#
DATA

#Examples:
#assign a red colored strip to leaf 9606, with label 'Human' (label is displayed in the mouseover popups)
#9606 #ff0000 Human

#assign a green, semi-transparent (alpha 0.5) strip to an internal node, without any label. If 'Show internal values' is set to 'No', this will only be displayed if the node is collapsed. 
#9606|5664 rgba(0,255,0,0.5)
";

}

sub popfontfile
{
print OFIL2 "TREE_COLORS
#use this template to define branch colors and styles, colored ranges and label colors/font styles
#lines starting with a hash are comments and ignored during parsing

#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file (except in the SEPARATOR line, which uses space).

#SEPARATOR TAB
#SEPARATOR SPACE
SEPARATOR COMMA

#First 3 fields define the node, type and color
#Possible types are:
#'range': defines a colored range (colored background for labels/clade)
#'clade': defines color/style for all branches in a clade
#'branch': defines color/style for a single branch
#'label': defines font color/style for the leaf label

#The following additional fields are required:
#for 'range', field 4 defines the colored range label (used in the legend)

#The following additional fields are optional:
#for 'label', field 4 defines the font style ('normal',''bold', 'italic' or 'bold-italic') and field 5 defines the numeric scale factor for the font size (eg. with value 2, font size for that label will be 2x the standard size)
#for 'clade' and 'branch', field 4 defines the branch style ('normal' or 'dashed') and field 5 defines the branch width scale factor (eg. with value 0.5, branch width for that clade will be 0.5 the standard width)

#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages
#=================================================================#
#        Actual data follows after the DATA keyword               #
#=================================================================#
DATA
#NODE_ID,TYPE,COLOR,LABEL_OR_STYLE,SIZE_FACTOR

#Examples
#internal node with solid branches colored blue and twice the standard width
#9031|9606,clade,#0000ff,normal,2
#internal node with dashed branches colored red and one half the standard width
#601|340,clade,#ff0000,dashed,0.5
#a single internal branch colored green, dashed and 5 times the normal width
#915|777,branch,#00ff00,dashed,5

#colored range covering all leaves of an internal node, colored red and with label 'Eukaryota'
#184922|9606,range,#ff0000,Eukaryota
#examples of colored ranges from iTOL's Tree of Life
#2190|2287,range,#aaffaa,Archaea
#623|1502,range,#aaaaff,Bacteria

#leaf label for node 9606 will be displayed in green, bold and twice the regular font size
#9606,label,#00ff00,bold,2

#leaf label for node 9031 will be displayed in yellow, bold italic and half the regular font size
#9031,label,#ffff00,bold-italic,0.5

#leaf label for node 8015 will be displayed in blue
#8015,label,#0000ff
";
}

sub popcolhash
{

#good general color scheme
%color_hash = (
	     "Actinobacteria","#2D33E7 Bacteria",
	     "Firmicutes","#2D33E7 Bacteria",
	     "Proteobacteria","#2D33E7 Bacteria",
		 "other_Bacteria","#2D33E7 Bacteria",
		 "Aquificae","#2D33E7 Bacteria",
		 "Armatimonadetes","#2D33E7 Bacteria",
		 "BacteroidetesChlorobi","#2D33E7 Bacteria",
		 "Caldiserica","#2D33E7 Bacteria",
		 "ChlamydiaeVerrucomicrobia","#2D33E7 Bacteria",
		 "Chloroflexi","#2D33E7 Bacteria",
		 "Chrysiogenetes","#2D33E7 Bacteria",
		 "Cyanobacteria","#2D33E7 Bacteria",
		 "Deferribacteres","#2D33E7 Bacteria",
		 "DeinococcusThermus","#2D33E7 Bacteria",
		 "Dictyoglomi","#2D33E7 Bacteria",
		 "Elusimicrobia","#2D33E7 Bacteria",
		 "FibrobacteresAcidobacteria","#2D33E7 Bacteria",
		 "Fusobacteria","#2D33E7 Bacteria",
		 "Gemmatimonadetes","#2D33E7 Bacteria",
		 "Nitrospinae","#2D33E7 Bacteria",
		 "Nitrospirae","#2D33E7 Bacteria",
		 "Planctomycetes","#2D33E7 Bacteria",
		 "Spirochaetes","#2D33E7 Bacteria",
		 "Synergistetes","#2D33E7 Bacteria",
		 "Tenericutes","#2D33E7 Bacteria",
		 "Thermodesulfobacteria","#2D33E7 Bacteria",
		 "Thermotogae","#2D33E7 Bacteria",
		 "other_Fungi","#28b1aa Fungi",
		 "Microsporidia","#28b1aa Fungi",
		 "other_Pezizomycotina","#28b1aa Fungi",
		 "Dothideomycetes","#28b1aa Fungi",
		 "Eurotiomycetes","#28b1aa Fungi",
		 "Sordariomycetes","#28b1aa Fungi",
		 "Leotiomycetes","#28b1aa Fungi",
		 "Taphrinomycotina","#28b1aa Fungi",
		 "Basidiomycota","#28b1aa Fungi",
		 "Mucoromycota", "#28b1aa Fungi",
		 "other_Dikarya","#28b1aa Fungi",
		 "other_Ascomycota","#28b1aa Fungi",
		 "Saccharomycotina","#e31a1c Saccharomycotina",
);

}