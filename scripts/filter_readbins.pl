#!/usr/bin/env perl

use strict;
use warnings;

if (@ARGV < 2)
{
	die "1.read bins file (use - for reading from STDIN) - 2.file with chromosomes to keep in first column  - 3.optional: ignore chr prefix and .fa suffix\n";
}

my $binsfile= shift;
my $chromfile = shift;
my $ignore = shift;

open (CF, $chromfile) or die "could not open chromlength file $chromfile: $!\n";
open (BF, $binsfile) or die "could not open read bins file $binsfile: $!\n";

my @fields = ();
my %chroms = ();
my $chr;
while (<CF>)
{
	if ($_ =~ /^#/)
	{
		next;
	}
	chomp;
	@fields = split ("\t", $_);
	# get rid of possible chr and .fa if wanted
	if ($ignore)
	{
		($chr = $fields[0]) =~ s/^chr//;
		$chr =~ s/^\.fa$//;
		$chroms{$chr} = 1;
	}
	else
	{
		$chroms{$fields[0]} = 1;
	}
}
close CF;
print STDERR "chromosomes to keep: ";
foreach $chr (sort keys %chroms)
{
	print STDERR "$chr ";
}
print STDERR "\n";

my $line = "";
my $all = 0;
my $kept = 0;
# read bins: chrom is in first column
while ($line=<BF>)
{
	$all++;
	@fields = split ("\t", $line);
	if ($ignore)
	{
		($chr = $fields[0]) =~ s/^chr//;
		$chr =~ s/^\.fa$//;
		if (exists $chroms{$chr})
		{
			print $line;
			$kept++;
		}
	}
	else
	{
		if (exists $chroms{$fields[0]})
		{
			print $line;
			$kept++;
		}
	}
}
close BF;
print STDERR "from $all lines, kept $kept with selected chromosomes\n";
if ($kept < 1)
{
	die "Error in filtering read bins (filter_readbins.pl): no lines were extracted from input. Probably the prefixes of the chromosome file and the BAM file are not compatible\n";
}
exit;
