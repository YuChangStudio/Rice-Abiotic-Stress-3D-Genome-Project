#!/public/home/ychang/anaconda3/bin/perl
use strict;
use Getopt::Long;
use List::Util 'sum';
use vars qw($list $gff $o);
Getopt::Long::GetOptions(
    'list=s'    => \$list,
    'o=s'   => \$o,
);

my %tpm;
my @samples;
my $count = 0;
open(IN, $list);
while(<IN>)
{
	chomp;
	my @arr = split(/ /,$_);
	$samples[$count] = $arr[0];
	my $tab = $arr[1];
	open(TAB,$tab);
    <TAB>;
	while(<TAB>)
	{
		chomp;
		my @arr = split(/\t/, $_);
        my $gid = $arr[0];
		$tpm{$gid} -> {$samples[$count]}=$arr[8];
	}
	close(TAB);
	$count+=1;
}
close(IN);

open(OUT,">$o");
print OUT "gene_id"."\t".join("\t", @samples)."\n";
foreach my $key (keys %tpm)
{
	my @arr;
	for(my $i = 0; $i <= $#samples; $i++)
	{
	 	if(defined $tpm{$key}->{$samples[$i]})
	 	{
	    	push(@arr, $tpm{$key}->{$samples[$i]});
	 	}
	 	else
	 	{
	   	 push(@arr, 0);
	 	}
	}
print OUT $key."\t".join("\t", @arr)."\n";
}
close(OUT);