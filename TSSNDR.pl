use strict;
use warnings;

my $f = shift;

my %dup;
my %h;
open F,"zcat $f|";
while(<F>){
	chomp;
	my @a = split;
	my ($gene,$nm,$region) = split /:/,$a[-2];
	next if ++$dup{$gene}> 3;
	$h{"$gene:$nm"}{$region} = $a[-1];
}
close F;

for my $gene (sort keys %h){
	my $nf_score = '0.00000';
	if( $h{$gene}{'TSS1'} + $h{$gene}{'TSS2'} > 0){
		$nf_score = sprintf "%.5f",$h{$gene}{'NDR'} / ($h{$gene}{'TSS1'} + $h{$gene}{'TSS2'}) / 2;
	}
	print "$gene\t$nf_score\n";
}
