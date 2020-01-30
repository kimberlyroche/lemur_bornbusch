use strict;
use warnings;
use List::MoreUtils qw(uniq);

my @tax_strings = qw();
open(DATA, "Bornbusch_16S_Y1_fecal_count-table.tsv");
while(<DATA>) {
	if($_ =~ /(D_0__.*?)\t/) {
		push @tax_strings, $1;
	}
	if($_ =~ /Unassigned/) {
		print($_);
	}
}
close DATA;

my @unique_strings = uniq @tax_strings;

print($#unique_strings."\n");
