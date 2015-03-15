#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Graph;

##Globals
my $scriptname		= $0;
my $VERSION			= '0.3';
my $parse			= q{};
my $metabolic		= q{};
my $regulatory		= q{};
my $ensembl			= q{};
my $buildm			= q{};
my $buildrpm		= q{};

##Handle arguments
if (@ARGV < 1){
	print "\n Try '$scriptname --man' for full info\n\n";
	exit(0);
}
else {
	GetOptions(	'help'			=> sub { pod2usage(1); },
				'version'		=> sub {print STDOUT "\n $scriptname version $VERSION\n"; exit(0) },
				'man'			=> sub {pod2usage(-exitstatus => 0, -verbose => 2); },
			  	'parse:i'		=> \$parse,
			  	'metabolic'		=> \$metabolic,
			  	'regulatory'	=> \$regulatory,
			  	'ensembl'		=> \$ensembl,
				'buildrpm'		=> \$buildrpm
			  );
}

if ($ensembl){
	my $ensembl_map = $ARGV[0];
	my $spot_file = $ARGV[1];
	my $FILE1;
	my $FILE2;
	open($FILE1, "<", $ensembl_map) or die "$0: failed to open ensembl map file $ensembl_map :$!\n";
	my $microarray_graph = Graph::Undirected->new;
	my $line;
	my $gene;
	my $spot;
	my %genehash;
	my $neighbor;
	my $spotid;
	my @neighbors;
	my $regulation;
	my $count;
	while(!eof($FILE1)){
		$line = readline($FILE1);
		chomp $line;
		($gene, $spot) = split(/\s+/, $line);
		$microarray_graph -> add_edge($gene, $spot);
	}
	my @vertices = $microarray_graph->vertices;
	#foreach my $vert (@vertices) {print "$vert\n";}
	close $FILE1;
	open($FILE2, "<", $spot_file) or die "$0: failed to open spot file $spot_file :$!\n";
	while(!eof($FILE2)){
		$line = readline($FILE2);
		($spotid, $regulation, $count) = split(/\t/, $line);
		chomp $spotid;
		chomp $regulation;
		$spotid =~ s/ //g;
		if($microarray_graph->has_vertex($spotid)){
			@neighbors = $microarray_graph->neighbors($spotid);
			foreach $neighbor (@neighbors){
				$genehash{$neighbor} = $regulation;
			}
		}
	}
	if ($buildrpm){
		my %spothash = %genehash;
		my $grn_edge_file = $ARGV[2];
		my $ppi_edge_file = $ARGV[3];
		my $met_edge_file = $ARGV[4];
		my $REGOUT = $ARGV[5];
		my $PPIOUT = $ARGV[6];
		my $METOUT = $ARGV[7];
		buildrpm(\%spothash, $grn_edge_file, $REGOUT, "grn");
		buildrpm(\%spothash, $ppi_edge_file, $PPIOUT, "ppi");
		buildrpm(\%spothash, $met_edge_file, $METOUT, "hmn");
	}
}

if ($metabolic){
	my $node_file = $ARGV[0];
	my $edge_file = $ARGV[1];
	my $graph = Graph::Undirected->new;
	my $lookup_graph = Graph::Undirected->new;
	my $FILEE;
	my $FILEN;
	my $node1;
	my $node2;
	my $line;
	open($FILEE, "<", $edge_file) or die "$0: failed to open edge file $edge_file :$!\n";
	while(!eof($FILEE)){
		$line = readline($FILEE);
		chomp $line;
		($node1, $node2) = split(/\s+/, $line);
		$graph -> add_edge($node1, $node2);
	}
	close($FILEE);
	my $rxn;
	my $genel;
	open($FILEN, "<", $node_file) or die "$0: failed to open node file $node_file :$!\n";
	while (!eof($FILEN)){
		$line = readline($FILEN);
		chomp $line;
		($rxn, $genel) = split(/\s+/, $line);
		$lookup_graph -> add_edge($rxn, $genel);
	}
	close($FILEN);
	my %existential_hash;
	my @vertices = $graph->vertices;
	my $ex1;
	my $ex2;
	foreach my $vertex (@vertices){
		my @lookup = $lookup_graph->neighbors($vertex);
		my @neighbors = $graph->neighbors($vertex);
		foreach my $neighbor (@neighbors){
			my @lookup_neighbor = $lookup_graph->neighbors($neighbor);
			foreach my $i (@lookup){
				foreach my $j (@lookup_neighbor){
					$ex1 = "$i\t$j\n";
					$ex2 = "$j\t$i\n";
					unless (($i eq $j) or exists($existential_hash{$ex1}) or exists($existential_hash{$ex2})){
						print $ex1;
						$existential_hash{$ex1} = 1;
						$existential_hash{$ex2} = 1;
					}
				}
			}
		}
	}
}

if ($regulatory){
	my $map_file = $ARGV[0];
	my $regulatory_network = $ARGV[1];
	my $lookup = Graph::Undirected->new;
	my $FILE;
	my $line;
	my $node1;
	my $node2;
	open($FILE, "<", $map_file) or die "$0: failed to open map file $map_file :$!\n";
	while(!eof($FILE)){
		$line = readline($FILE);
		chomp $line;
		($node1, $node2) = split(/\t/, $line);
		$lookup -> add_edge($node1, $node2);
	}	
	close($FILE);
	open($FILE, "<", $regulatory_network) or die "$0: failed to open map file $regulatory_network :$!\n";
	while(!eof($FILE)){
		$line = readline($FILE);
		chomp $line;
		($node1, $node2) = split(/\t/, $line);
		my @neighbors1 = $lookup->neighbors($node1);
		my @neighbors2 = $lookup->neighbors($node2);
		foreach my $neighbor1 (@neighbors1){
			foreach my $neighbor2 (@neighbors2){
				print "$neighbor1\t$neighbor2\n";
			}
		}
	}	
	close($FILE);
	
}

elsif ($parse){
	my $expression_file = $ARGV[0];
	my $network_file = $ARGV[1];
	my $expression_hash_ref = "";
	$expression_hash_ref = parse_expression_file($expression_file, $parse);
	create_subnetwork($network_file, $expression_hash_ref);
}

sub parse_expression_file{
	my ($filename, $parse_number) = @_;
	my %expression_hash = ();
	my $expression = "";
	my $line = "";
	my @line_array = ();
	my $FILE;
	open($FILE, "<", $filename) or die "$0 : failed to open expression file $filename : $!\n";
	while(!eof($FILE)){
		$line = readline($FILE);
		chomp $line;
		@line_array = split(/\s+/, $line);
		$expression = $line_array[$parse_number];
		$expression_hash{$expression} = 1;
	}
	close($FILE);
	return \%expression_hash;
}

sub create_subnetwork{
	my $filename = shift;
	my $expression_hash_ref = shift;
	open(my $FILE, "<", $filename) or die "$0 : failed to open edge file $filename : $!\n";
	my $line = "";
	my @line_array = "";
	my $edge1 = "";
	my $edge2 = "";
	my %expression_hash = %$expression_hash_ref;
	while(!eof($FILE)){
		$line = readline($FILE);
		chomp $line;
		($edge1, $edge2) = split(/\s+/, $line);
		if(exists $expression_hash{$edge1} and exists $expression_hash{$edge2}){
			print $line."\n";
		}
	}
	close($FILE);
}

sub buildrpm{
	my $expression_hash_ref = shift;
	my $infile = shift;
	my $outfile = shift;
	my $netname = shift;
	my $node1;
	my $node2;
	my $line;
	my %expression_hash = %$expression_hash_ref;
	open(my $FILEI, "<", $infile) or die "$0 : failed to open file $infile :$!\n";
	open(my $FILEO, ">", $outfile);
	while(!eof($FILEI)){
		$line = readline($FILEI);
		chomp $line;
		($node1, $node2) = split(/\s+/, $line);
		if ($node1 ne $node2){
			if (exists $expression_hash{$node1} and exists $expression_hash{$node2}){
				print $FILEO $node1."\t".$node2."\t".$expression_hash{$node1}."\t".$expression_hash{$node2}."\n";
			}
		}
	}
	close $FILEI;
	close $FILEO;
}
## POD documentation
=pod

=head1 NAME

integrate_subnetworks.pl

=head1 VERSION

Documentation for integrate_subnetworks.pl version 1.0

=head1 SYNOPSIS

integrate_subnetworks.pl [--parse=<NUMBER>] EXPRESSION_FILE EDGE_FILE [>OUTPUT]
integrate_subnetworks.pl [--metabolic]  NODE_FILE EDGE_FILE [>OUTPUT]
integrate_subnetworks.pl [--regulatory] MAP_FILE REGULATORY_FILE [>OUTPUT]
integrate_subnetworks.pl [--ensembl]	ENSEMBL_MAP EXPRESSION_FILE [--buildrpm] REGULATORY_NETWORK PPI_NETWORK MET_NETWORK REG_OUT PPI_OUT MET_OUT

=head1 DESCRIPTION

Script for generating subnetworks from expression files and a given network.

=head1 OPTIONS

Mandatory arguments

=over 8

=item B<-e, --ensembl>
Parses a ensembl to affymetrix map.

=item B<-p, --parse=>I<number>
In the expression file, extract the nth item.

=item B<-m, --metabolic>
Parses metabolic file.

=item B<-r, --regulatory>
Parses regulatory file.

=item B<-h, --help>
Prints help messgae and exits.

=item B<--man>
Displays the manual page.

=item B<-o, --outfile=>I<filename>
Prints directly to file I<filename> instead of standard out.

=item B<-buildrpm --r>R<FILENAME> P<FILENAME> M<FILENAME> R_OUT<FILENAME> P_OUT<FILENAME> M_OUT<FILENAME>
Builds subnetworks for regulator, protein interaction and metabolic interactions, based on the existing expression.

=back

=head1 USAGE

Examples:
	
	integrate_subnetworks.pl --parse=2 gastric_ensembl.list edgefile_unique_EnsEMBL50_PPI_pairs_w_blasthits_wo_selfints.list > gastric_ppi_subnetwork.list
	integrate_subnetworks.pl -p=2 gastric_ensembl.list edgefile_unique_EnsEMBL50_PPI_pairs_w_blasthits_wo_selfints.list > gastric_ppi_subnetwork.list
	integrate_subnetworks.pl --metabolic palsson_nodes.txt palsson_edges_50_new.txt
	integrate_subnetworks.pl --regulatory ensembl_uniprot_map.txt human_gene_regulatory_network.txt
	integrate_subnetworks.pl --ensembl ensembl_affyHG_U133_PLUS_2.txt
	integrate_subnetworks.pl --ensembl ensembl_affyHG_U133_PLUS_2.txt spots_in_prostrate.txt  --buildrpm grn_edges.list ppi_edges.list gene_metabolic_network_10.txt grn_prostrate.out ppi_prostrate.out hmn_prostrate.out


=head1 AUTHOR

Written by Corey M. Hudson

=head1 DEPENDENCIES

Uses Perl Modules Getopt::Long and Pod::Usage and Graph

=head1 LICENSE AND COPYRIGHT

Copyright (c) 2012 Corey M. Hudson.
All rights reserved.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details. 
http://www.gnu.org/copyleft/gpl.html 


=cut
