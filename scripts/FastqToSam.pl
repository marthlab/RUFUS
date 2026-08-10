#!/usr/bin/perl
# Convert a FASTQ stream to unmapped SAM records on stdout, for RUFUS's whole-genome FASTQ input path.
# Usage:  FastqToSam.pl <fastq> [header]
#
# Used ONLY to feed k-mer COUNTING: the RUFUS generator pipes this through `samtools fastq` into
# jellyfish, which just needs the sequences — read pairing is irrelevant here, so every record is
# emitted unmapped (flag 4, MAPQ 0). Paired filtering is done separately from the raw FASTQs via the
# -q1/-q2 path; single-end uses RUFUS.Filter.single. The previous MAPQ='*' + trailing tab made samtools
# reject the stream. The generator emits ONE @HD header for the whole stream (samtools rejects a header
# mid-stream), so only the first FastqToSam.pl call in a generator is passed 'header'.
use strict;
use warnings;

my $emit_header = (defined $ARGV[1] && $ARGV[1] eq 'header');

open(my $Fastq, '<', $ARGV[0]) || die "ERROR could not open fastq file $ARGV[0]";

print "\@HD\tVN:1.6\tSO:unsorted\n" if $emit_header;

while (my $l1 = <$Fastq>) {
    my $l2 = <$Fastq>;   # sequence
    my $l3 = <$Fastq>;   # '+'
    my $l4 = <$Fastq>;   # quality
    last unless defined $l4;
    chomp $l1;
    chomp $l2;
    chomp $l4;
    my @t   = split ' ', $l1;
    my $name = substr($t[0], 1);   # strip leading '@'
    # QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL   (unmapped: flag 4, MAPQ 0)
    print "$name\t4\t*\t0\t0\t*\t*\t0\t0\t$l2\t$l4\n";
}
close($Fastq);
