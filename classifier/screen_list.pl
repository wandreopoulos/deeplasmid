#!/usr/bin/env perl
use warnings;


$n_args = @ARGV;
if (($n_args != 2) && ($n_args != 3)) {
    print "Filter a fasta file for/against a list of entry names.\n";
    print "Usage: ./screen_list.pl <list> <fasta file> <<keep?>>\n";
    print "If a third argument is given the list entries will be kept,\n";
    print "otherwise they will be rejected.\n";
	print "Also, only unique IDs will be printed out...duplicates discarded.\n";
    exit;
}

my $invert_sense = 0;
if ($n_args == 3) {$invert_sense = 1;}

#
# open the list file
#
open(F,$ARGV[0]) || die "Couldn't open file $ARGV[0]\n";
my %bad_ones = (); 
while (my $i = <F>) {
    chomp $i;
    $i =~ s/\s*$//g;
    $bad_ones{$i} = 1;
}
close F;

#
# open fasta file
#
open(F,"cat $ARGV[1] |") || die "Couldn't open file $ARGV[1]\n";
my %seen=();
my $good_one = 1;
my $id;
my $print;
my ($counter,$duplicates)=(0,0);
while (my $i = <F>) {
    if ($i =~ />\s*(\S+)/) {
	    $id = $1;
	    if (exists($bad_ones{$id})) {
            $good_one = 0;
			delete $bad_ones{$id};
			$counter++;
        }else {
            $good_one = 1;
        }

		if (exists $seen{$id}){
			$print='no';
			$duplicates++;
		}else{
			$seen{$id}='';
			$print='yes';
		}
    }

	next if $print eq 'no';

    if (($invert_sense == 0) &&  ($good_one == 1)) {print $i;}
    elsif (($invert_sense == 1) &&  ($good_one == 0)) {print $i;}
}
close F;

#foreach (keys %bad_ones){
#	print STDERR "left over: $_\n";
#}
#print STDERR "sequences printed: $counter\n";
#print STDERR "sequences duplicated: $duplicates\n"; # this is wrong
#print "Line Count:\t$line_counter\n";
#print "Read Count:\t$read_counter\n";
