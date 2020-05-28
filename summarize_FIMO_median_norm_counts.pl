use strict;
use warnings;
use Statistics::PointEstimation;
use Statistics::TTest;
use Statistics::Basic qw(:all);

#sub average
#{
#	my @array = @_;
#	my $sum;
#	foreach (@array)
#	{
#		$sum += $_;
#	}
#	return $sum/@array;
#}
if (scalar(@ARGV) != 2)	
{
	die ("Usage: $0 <Input file with per sample normalized FIMO counts> <Output file with summarized average FIMO counts>");
}
	my %merge_normal=();
	my %merge_tumor=();

	my $info = $ARGV[0];
	my @info_array = split(/\./,$info);
	my $gene = $info_array[3];
	my $miR = $info_array[4];
	
	open (FH, "$ARGV[0]");
	open (OUT, ">$ARGV[1]");
#	print OUT "target_gene\tconserved_miR\tNormal_normalized_average\tTumor_normalized_average\n";
	while(my $l = <FH>)	
{
		chomp $l;
		my @a = split(/\t/,$l);
		if ($a[1] =~ m/Normal/)
		{
			push (@{$merge_normal{$gene}},$a[3]);
		}
		elsif ($a[1] =~ m/Tumor/)
		{
			push (@{$merge_tumor{$gene}},$a[3]);
		}
	}
	close FH;
foreach my $key (sort keys %merge_normal)
{
	my @normal=@{$merge_normal{$key}};
	my @tumor=@{$merge_tumor{$key}};
	my $ttest = new Statistics::TTest;
#	$ttest->set_significance(90);
	$ttest->set_significance(95);
	$ttest->load_data(\@normal,\@tumor);
	my $t = $ttest->t_statistic;
	my $prob = $ttest->{t_prob};
	my $test = $ttest->null_hypothesis();
#	my $average_normal = sprintf("%.3f", average(@normal));
#	my $average_tumor = sprintf("%.3f", average(@tumor));
	my $median_normal = sprintf("%.3f", median(@normal));
	my $median_tumor = sprintf("%.3f", median(@tumor));
	my $out_t = sprintf("%.3f", $t);
	my $out_prob = sprintf("%.3f", $prob);
	if ($median_normal == 0 && $median_tumor == 0)
	{
		 print OUT "$gene\t$miR\t$median_normal\t$median_tumor\t-\t-\t-\n";
	}
	else
	{
		print OUT "$gene\t$miR\t$median_normal\t$median_tumor\t$out_t\t$out_prob\tNull hypothesis is $test\n";
	}
}

	close OUT;




