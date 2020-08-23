#!/usr/bin perl
use strict;
use warnings;

use JSON;
use Data::Dumper;
use JSON::Create 'create_json';
use List::MoreUtils qw(uniq);

open FILE, $ARGV[0] or die "file does not exsist!\n";

#pmid	section	sentence_number	geneId	geneoffset	diseaseId	diseaseoffsets	sentence

my (%sentence, $sentence, %gene, $gene, %disease, $disease);
my $gene_tag = "";
my $disease_tag = ""; 

while (my $_ = <FILE>) {
	chomp $_;
	my @array = split (/\t/, $_);
	my $array;
	$array[0] =~ s/^ //g;
	$array[2] =~ s/^ //g;
	my $key = $array[0]."-".$array[2];
	push (@{$sentence{$key}}, $array[7]);
	my $gene_flag = $array[4]."#gene".$array[3];
	if ($gene_flag ne $gene_tag) {
		push (@{$gene{$key}}, $gene_flag); 
		$gene_tag = $gene_flag;
	}	
	my $disease_flag = $array[6]."#disease".$array[5];
	if ($disease_flag ne $disease_tag) {
		push (@{$disease{$key}}, $disease_flag);
		$disease_tag = $disease_flag;
	}
}

for $sentence  ( keys %sentence) {
	my ($subj, $obj);
	my @subj = &split($sentence, @{$gene{$sentence}});
	my @obj = &split($sentence, @{$disease{$sentence}});
	my @denote = (@subj, @obj);
	my @relation = &combi($sentence, \@{$gene{$sentence}}, \@{$disease{$sentence}});
	my $pmid;
	my @pmid = split (/-/, $sentence);
	my %example = (
		text => uniq(@{$sentence{$sentence}}),
		sourcedb => "PubMed",
		sourceid => $pmid[0],
		denotations => \@denote,
		relations => \@relation,
	);
	open OUT, ">$ARGV[1]/$sentence.json" or die;
	my $json = create_json (\%example);
	print OUT $json;
	close OUT;
#	print Dumper \%example;
#	my $json = create_json (\%example);
#	print $json;
}



sub split() {
	my ($count, @text) = @_;
	my @out;
	foreach my $text (@text) {
		my $array;
		my @array = split (/#/, $text);
		my %out = (
			id => $count."#".$text,
			span => {
				begin => $array[0],
				end => $array[1],
			},
			obj => $array[2],	
		);
		push @out, {%out};
	}
	return @out;
}

sub combi() {
	my ($_, @arr) =@_;
	my @subj = @{$_[1]};
	my @obj = @{$_[2]};
	my @out;
	for (my $i=0; $i<= $#subj; $i++)  {
		for (my $j=0; $j<= $#obj; $j++) {
			my %out = (
				pred => "associated_with",
				subj => $_ . "#".$subj[$i],
				obj => $_ . "#" .$obj[$j],
				id => $subj[$i].$obj[$j],
			);
			push @out, {%out};
		}
	}	
	return @out;

}
