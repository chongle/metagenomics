#!/usr/bin/perl -w

# converts a multi line sequence in a fasta file into one line sequence... easier to use
# output is printed on screen, so use pipe to save in a file
# also the sequence is converted to upper case

if(@ARGV!=1) { die "correct usage:\n\t perl format_fasta.pl <fasta_file> "};

open(FASTA,$ARGV[0]) || die "can't open fasta file";
$line="";
while($line !~ /^>/){
	$line=<FASTA>;
}
$line  =~ s/[\r]//g; # remove carriage return
print $line;
$prev="";
while (1){
	$line=<FASTA>;
	$line=~ s/[\r]//g; # remove carriage return
	while ($line !~ /^>/){
		chomp $line;
		$prev=uc $prev.$line;
		$prev =~ s/[\r]//g; # remove carriage return
		$line=<FASTA>;
		if(!(defined $line)) { print $prev."\n"; close(FASTA); exit(1);}
	}
	print $prev."\n";
	$line  =~ s/[\r]//g; # remove carriage return
	print $line;
	$prev="";
}


