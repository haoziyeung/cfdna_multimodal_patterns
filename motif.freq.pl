use strict;
use warnings;
#use Data::Dumper;
use List::Util 'sum';
use feature 'say';
#$Data::Dumper::Sortkeys = 1;


my ($infile,$featureFile) = @ARGV;
die "perl $0 data/PCSC513N2Y19KB096C.gz abnormalHap.final.list\n\n" unless $infile && -e $infile && $featureFile && -e $featureFile;


my %h;
my $total;
open F,"zcat $infile|";
while(<F>){
    chomp;
    my @a = split;
	$h{$_} ++;
	$total ++;
}
close F;

my @features;
open F,$featureFile;
while(<F>){
    chomp;
	my $freq = '0.0000000000';
	if(exists $h{$_}){
		$freq = sprintf "%.10f",$h{$_} /  $total;
	}
	print "$_\t$freq\n";
}
close F;


