#!/usr/bin/perl
use strict;
use warnings;
my @files = glob("*.tree");
 foreach my $file ( @files ) {
  next unless $file=~ m/(?<ID>.*).tree/gi;
  my $id_name=$+{ID};
my $out = "itol_run.bat";
open (OUT,">>", $out) or die "can't write the file $out:$!";
print OUT "perl create_iTOL_tree.pl $file Yeast_HGT\n";
close OUT;
}