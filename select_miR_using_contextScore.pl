
use strict;
use warnings;

 if (scalar(@ARGV) != 2)	
 {
 die ( "Usage:\nselect_miR_using_contextScore.pl <Input file containing site type number> <Ouput file>\n" );
 }

my %hash1=();

#Form a unique key with chr_pos from merged list.report
open (FIRST,"$ARGV[0]");
open (OUT1, ">$ARGV[1]");
my $firstLine = <FIRST>;
print OUT1 $firstLine;
while(my $line = <FIRST>)	
{
	next if ($line =~ m/seed_type/);
	chomp $line;
	my @array = split('\t',$line);
	my $uniq = join("\t",$array[0],$array[2],$array[3],$array[4],$array[5]);
	my $value = join ("\t",$array[1],$array[6]);
	push( @{$hash1{$uniq}},$value);
}
close FIRST;

foreach my $key ( keys %hash1)
{
	my @all = @{$hash1{$key}};
	my @keys = split ("\t",$key);
	my $count = scalar(@all);
	if ($count == 1)
	{
		my @final = split ("\t",$all[0]);
		print OUT1 "$keys[0]\t$final[0]\t$keys[1]\t$keys[2]\t$keys[3]\t$keys[4]\t$final[1]\n";
	}
	elsif ($count > 1)
	{
		my @final = split ("\t",$all[0]);
		print OUT1 "$keys[0]\t$final[0]\t$keys[1]\t$keys[2]\t$keys[3]\t$keys[4]\t$final[1]\n";
	}
}	
close FIRST;
close OUT1;

