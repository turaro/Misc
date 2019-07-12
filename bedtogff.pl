#Parse pseudogff to add exons to corresponding genes
use strict;
use warnings;

my $i;
my @array;
my $type;
my $gene_id;
my $mrna_id;
my $exon_id;
my $cds_id;
my $pure_id;
open GFF, "<$ARGV[0]" or die $!;
while ( my $line = <GFF> ) {
	chomp $line;
    @array = split( "\t", $line );
    $type = $array[2];
		if ($type eq "gene") {
			$array[3]+=1;
			$gene_id = $array[8];
			$i = 1;
			print join "\t", @array, "\n";
			$array[2]="mRNA";
			$array[8] =~ /ID=(.*)/;
			$pure_id = $1;
			$array[8] = join "", "ID=", $1, ".1;Parent=", $1;
			$mrna_id = $array[8];
			print join "\t", @array, "\n";
		}
		elsif ($type eq "exon") {
			$array[3]+=1;
			$exon_id = join ".1:exon:", $gene_id, $i;
			$exon_id = join "", $exon_id, ";Parent=", $pure_id, ".1";
			push @array, $exon_id;
			print join "\t", @array, "\n";
			pop @array;
			$array[2]="CDS";
			$cds_id = join ".1:cds:", $gene_id, $i;
			$cds_id = join "", $cds_id, ";Parent=", $pure_id, ".1";
			push @array, $cds_id;
			print join "\t", @array, "\n";
			$i++;
		}
	
	}
