#!/usr/bin/perl

# extended flagstats with number of secondary and supplementary alignments,
# paired end reads mapping to different chromosomes,
# bucketsort insert sizes, calculate average, median, standard deviation and "SDpercent" as standard deviation/median*100
# print out the insert size frequencies on stdout for R:
# insertsize_bin\tnumber_reads
# optional, print out median and "SDpercent" to a file (mimicking the text output of insertSizeDistribution.r)

use strict;
use warnings;
use Getopt::Std;

my %opts = (n=>1);
getopts("c:i:o:n", \%opts);

my $samfile = $opts{i};
my $peaberrstat = $opts{o};
my $chromfile = $opts{c};
my $minmapq = $opts{n};

if (defined $opts{h} || ! defined $samfile || ! defined $chromfile || ! defined $peaberrstat)
{
	die "USAGE: $0 [options]
	-c FILE file with the chromosomes that are wanted (must fit w.r.t. pre- and suffix!)
	-i FILE	SAM file (set to - for pipe from STDIN)
	-n INT minimal mapping quality for improper pair (default 1)
	-o FILE file to write the percentage of improper pairs out to (_DiffChroms.png_qcValues.txt)
	-h help\n";
}

open (IN, $samfile) or die "Cannot open $samfile: $!\n";
open (CF, $chromfile) or die "Cannot open $chromfile: $!\n";
open (PE, ">$peaberrstat") or die "could not open $peaberrstat for writing: $!\n";

if ($minmapq !~ /^\d+$/)
{
	print STDERR "minimal mapping quality $minmapq is not a number, setting to default 1\n";
	$minmapq = 1;
}

# first read in which chromosomes should be regarded for the statistics
my @help = ();
my %chroms = ();	# to look up later if the chromosome is wanted
my @chromarray = ();	# for the output in matrix form
my $chrnum = 0;
while (<CF>)
{
	if ($_ =~ /^#/)
	{
		next;
	}
	chomp;
	@help = split ("\t", $_);
	$chroms{$help[0]} = 1;
	$chromarray[$chrnum] = $help[0];
	$chrnum++;
}
close CF;
print STDERR "$chrnum chromosomes to keep: @chromarray\n";

# extended flagstats
# all reads
my $all = 0;
# only the non-dup, non-secondary, non-supplementary reads
my $uniq = 0;
# of these, with min. mapping quality
my $minmapuniq = 0;
# of these, those on wanted chroms
my $onchr = 0;
# of these, both on wanted chroms - these will be the denominator
my $both = 0;
# the PE aberrations among them
my $aberrant = 0;
my %chrompairs = ();	# hash of hashes: chromosomes of read => chr. of mate => number
my $flag = 0;


while (<IN>)
{
	if ($_ =~ /^\@/)	# there might be a SAM header
	{next;}
	$all++;
	@help = split ("\t", $_);
	$flag = $help[1];
	# not unmapped, no duplicate and no secondary/supplementary alignment, and mapqual >= X
	if (!($flag & 4) && !($flag & 1024) && !($flag & 256) && !($flag & 2048))
	{
		$uniq++;
		if ($help[4] >= $minmapq)
		{
			# of these, with mapqual >= X
			$minmapuniq++;
		}
		else
		{
			next;
		}
		# is the read itself on a wanted chromosome?
		if (defined $chroms{$help[2]})
		{
			$onchr++;
			# and the mate on a wanted chromosome
			if (defined $chroms{$help[6]} || $help[6] eq "=")	# same chrom is usually indicated by "=" instead of repeating the name
			{
				$both++;
				# paired end aberration:  mate also has to be mapped, on a different chrom
				if (!($flag & 8) && $help[6] ne "=" && ($help[2] ne $help[6]))
				# keep matrix symmetrical to see whether there is a bias, e.g. more 1->10 than 10->1
				{
					$aberrant++;
					# only use read1 info, since read2 might have mapq 0, and the info of having bias w.r.t. which read
					# is more interesting
					if ($flag & 64)
					{
						$chrompairs{$help[2]}{$help[6]}++;
					}
				}
			}
		}
	}
}
close IN;

# in case we have single end reads, there will be no counted ones for isizes ($ctr < 1)
# and PE aberrations ($aberrant < 1) => fill the files with placeholder "NA"

if ($aberrant < 1)
{
	print STDERR "no aberrant paired reads found, single end reads?\n";
	print PE "NA\n";
}
else
{
	my $percentage = sprintf ("%2.2f\n", ($aberrant/$both*100));
	print PE $percentage;
}
close PE;

