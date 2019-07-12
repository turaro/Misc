#Script for extracting sequences by the list

print "Fasta file:\n";
my $name1=<STDIN>;
chomp $name1;
print "\n";

print "List file:\n";
my $name2=<STDIN>;
chomp $name2;
print "\n";

print "Output file:\n";
my $name3=<STDIN>;
chomp $name3;
print "\n";

open FILE, "<$name2" or die ("Can't open list file! :c $!\n\n"); 
my @words;
my @list;

while (my $line=<FILE>)  
 {
 chomp($line);
$line =~ s/^\s+//;
@words = split /\s+/, $line;
$words[0] = substr $words[0], 1, -1;
push @list, $words[0];
 }
close FILE;

my %storage;
my $content;

open(my $fh, '<', $name1) or die "Can't open fasta file! :c $!\n\n";
{
local $/;
$content = <$fh>;
}
close($fh);

my @fasta = split />/, $content;
my @two;

foreach my $i(@fasta) 
{
@two = split /\n/, $i, 2;
$two[0] =~ /.*gene=(.*)/;
$storage{$1} = $two[1];
}

my @output;

foreach my $j(@list) 
{
if (exists $storage{$j})
	{
open($output, ">>", "$name3") or die;
print $output ">$j\n$storage{$j}";
	}
}

