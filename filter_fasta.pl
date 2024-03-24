#!/usr/bin/perl
use strict;
use warnings;
my ($infile, $minlen) = @ARGV;

open my $IN, '<', $infile or die "Cannot open $infile: $!";

{
local $/=">";
while(<$IN>) {
chomp;
next unless /\w/;
my @keep = split /\n/;
my $header = shift @keep;
my $seqlen = length join "", @keep;
if($seqlen >= $minlen){
print ">$header\n", join("\n", @keep), "\n";
}
}
local $/="\n";
}

close $IN;
exit;


##script from https://bioinformaticsreview.com/20200728/extract-fasta-sequences-based-on-sequence-length-using-perl/
##corrected with chatgpt