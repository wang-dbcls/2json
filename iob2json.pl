#!/usr/bin/perl -w
require 5.000;

if ($#ARGV != 2) {die "\nperl *.pl iob seiob protein\nperl *.pl iob seiob sub\n\n";}

open INPUT, "<$ARGV[0]" or die;
open OUTPUT, ">$ARGV[1]" or die;
my $option= $ARGV[2];

while ($line=<INPUT>) {
	chomp $line;
	if ($line) {
		@info= split("	", $line);
		$timer=$#info-1;
		$word=$info[0];
		for ($i=1;$i<=$timer;$i++) {
			$word =join ("\t", $word, $info[$i]);
		}
		push (@word, $word);
		push (@tag, $info[-1]); #for aimed and genia tagged by the tagger trained with aimed
# for original genia
=pod
		if ($info[-1]=~ /DNA_domain_or_region/) {push (@tag, $info[-1]);}
		else {push (@tag, "O");}
=cut
	}	
	else {push (@word, '');push (@tag, '');}
}

if ($option eq '') {
	die 
	"\nPlease input option!\nprotein or sub\n\n";
}	
	

if ($option eq 'protein') {
	for ($i=0; $i<=$#tag; $i++) {
		$tag_tmp[$i] = $tag[$i];
		if (substr($tag[$i], 0, 1) eq "I") { #elsif A
			if (($i==$#tag)||(substr($tag[$i+1], 0, 1) ne "I")) {
				$tag_tmp[$i]= "E-protein";
			}
			else {
				$tag_tmp[$i]= "I-protein";
			}
		} #elsif A
		elsif (substr($tag[$i], 0, 1) eq "B") { #elsif B
			if (($i==$#tag)||(substr($tag[$i+1], 0, 1) ne "I")) {
				$tag_tmp[$i]= "S-protein";
			}
			else {
				$tag_tmp[$i]= "B-protein";
			}
		} #elsif B
		elsif (substr($tag[$i], 0, 1) eq "O") {
			substr($tag_tmp[$i], 0, 1) = "O";
		}
	
		if ($tag[$i] eq ''){
			print OUTPUT "\n";
		}
		else {
			print OUTPUT "$word[$i]	$tag_tmp[$i]\n";
		}
	}
}

if ($option eq 'sub') {
	for ($i=0; $i<=$#tag; $i++) {
		$tag_tmp[$i] = $tag[$i];
		if (substr($tag[$i], 0, 1) eq "I") { #elsif A
			if (($i==$#tag)||(substr($tag[$i+1], 0, 1) ne "I")) {
				substr($tag_tmp[$i], 0, 1) = "E"; #end of protein
			}
			else {
				substr($tag_tmp[$i], 0, 1) = "I";
			}
		} #elsif A
		elsif (substr($tag[$i], 0, 1) eq "B") { #elsif B
			if (($i==$#tag)||(substr($tag[$i+1], 0, 1) ne "I")) {
				substr($tag_tmp[$i], 0, 1) = "S"; #single protein word
			}
			else {
				substr($tag_tmp[$i], 0, 1) = "B";
			}
		} #elsif B
		elsif (substr($tag[$i], 0, 1) eq "O") {
			substr($tag_tmp[$i], 0, 1) = "O";
		}
	
		if ($tag[$i] eq ''){
			print OUTPUT "\n";
		}
		else {
			print OUTPUT "$word[$i]	$tag_tmp[$i]\n";
		}
	}
}

close INPUT;
close OUTPUT;
#!/usr/bin/perl

# Modified: Oct/01/2012
# one-line-output

open FILE, $ARGV[0] or die "perl *.pl seiob GENE\n";
open OUT, ">$ARGV[1]" or die "can't create result file [$ARGV[1]].";

while ($line=<FILE>) {
	chomp($line);	
	if ($line) {
		@info=split("\t", $line);
#		if ($info[0]=~ /\d{7,12}/) {next;}
		print OUT "$info[0]"."\|"."$info[-1] ";
	}
#	else {print OUT "\n";}
	else { print OUT " " }
}
close OUT;
close FILE;
#!/usr/bin/perl
require 5.000;
use warnings;
use strict;

# Nov/2016 (pir overlap annotation)
# Jan/2016 (subclass option removed, medline id 2 pmid [jnlpba])
# Oct/2015 (json output for pubannotation)
# Oct/2012 (offset counted on abstract level, including whitespace)

if ($#ARGV == -1) {die "\n$0 single-line-GENE raw-txt off-a1 json corpus\n\n"}

open FILE, "<$ARGV[0]" or die "can't open GENE file $!";
open RAW, ">$ARGV[1]" or die;
open OFF, ">$ARGV[2]" or die;
open JSON, ">$ARGV[3]" or die;
my ($text, $denotations, $class) = "";

# pmid
my ($filename, %hash) = "";
my @filename = split (/\//, $ARGV[0]);
my $each = substr ($filename[-1],0, -5);
my ($sourceid) = $filename[-1] =~ /(\d+)/;
my $id = $sourceid;
my $corpus = $ARGV[4];
if ( $corpus eq "jnlpba") {
	open LIST, "</home/wang/pubann/UID-PMID-PMCID.lst" or die;
	while (my $_ = <LIST>) {
		chomp $_;
		my @tab = split(/\t/, $_);
		my $tab;
		$hash{$tab[0]} = $tab[1];
	}
	close LIST;
	while (my ($key, $value) = each %hash ) {
		if ( $key = $id ) { $sourceid = $hash{$id} }
		else { print "Error, $id does not have a corresponding pmid.\n" }
	}
}


my $i = 0;
if ( $corpus eq "pir2") {
	open PIR, "/home/wang/pubannotation/data/pir/corpus2/protein/$sourceid.json" or die;
	my $pir = <PIR>;
	chomp $pir;
	my @c = $pir =~ /T\d+/g;
	$i = @c;
	close PIR;
}

my $separator = quotemeta ("|");	
while (my $fileline=<FILE>) {
	chomp($fileline);
	my $length = 0;
	my $start = 0;
	my $ending = 0;
	my $temp = 0;
	$i = $i + 1;
	my @info=split(" ", $fileline);	
	my $string = '';
#	print RAW "$corpus"."$. ";
	foreach my $info (@info) {
		my @word= split ($separator, $info);
		print RAW "$word[0] ";
		my @tag = split (/-/, $word[-1]);
		my $tag;
		$text .= &escape($word[0]) ." ";
		$length+= length($word[0])+1;
		if ($word[-1] =~ "-") { $class = substr($word[-1],2) }
		if ($word[-1] =~ "S-") {
			$start= $temp;
			$ending= $start+ length ($word[0]);
			if ($tag[-1] eq "gene") {
				print OFF "\tT$i $tag[-1] $start $ending\n";
				$denotations .= "{\"id\": \"T$i-$each-$start-$ending\", \"span\": \{\"begin\": $start, \"end\": $ending\}, \"obj\": \"$class\"},";
			}
			$i++;
		}
		elsif ($word[-1] =~ "B-") {
			$start=$temp;
			$string.=" $word[0]";
		}
		elsif ($word[-1] =~ "I-") {
			$string.=" $word[0]";
		}
		elsif ($word[-1] =~ "E-") {
			$ending=$temp+ length ($word[0]);
			$string.=" $word[0]";
			$string=~ s/^\s//;
			$string=~ s/\s$//;
			if ($tag[-1] eq "gene") {
				print OFF "\tT$i $tag[-1] $start $ending\n";
				$denotations .= "{\"id\": \"T$i-$each-$start-$ending\", \"span\": \{\"begin\": $start, \"end\": $ending\}, \"obj\": \"$class\"},";
			}
			$string="";
			$i++;
		}
		$temp=$length;
	}
	print RAW "\n";
}

if ( not defined $denotations) { $denotations = "," }
chop ($denotations);  
print JSON "{\"text\": \"$text\", \"sourcedb\": \"PubMed\", \"sourceid\": \"$sourceid\", \"denotations\": [$denotations]}";
close FILE;
close RAW;
close OFF;
close JSON;

sub escape {
	my $_ = shift;
	$_ =~ s/"/\\"/g;
	return $_;
}

