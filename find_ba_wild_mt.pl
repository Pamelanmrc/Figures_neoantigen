use strict;

# perl find_ba_wild_mt.pl  14-134_wild_Type_netmhcpan.output.txt    14-134_mutated_Type_netmhcpan.output.txt 14-134_Mutect2_somatic.PASS.maf VAF_file

#NBPF1_16565761_	HLA-A*02:01	VLQDSLDRC	0.0976360	1.919	0.283398	6.284	WB	NBPF1_16565761_	HLA-A*02:01	LLEAVEPEV	0.6734670	0.211	0.654098	0.573	SB

open(WRITE, ">SAMPLE_OUT");


`egrep \"<= SB|<= WB\" $ARGV[0] |  awk -F \" \" \'{ print \$11\"\t\"\$2\"\t\"\$10\"\t\"\$12\"\t\"\$13\"\t\"\$14\"\t\"\$15\"\t"\$16\"\t"\$18}\' > wild_BA`;
`grep \"<= SB\" $ARGV[1]  | awk -F \" \" \'{ print \$11\"\t\"\$2\"\t\"\$10\"\t\"\$12\"\t\"\$13\"\t\"\$14\"\t\"\$15\"\t\"\$16}\'  > mutant_BA`;


open(FILE, "mutant_BA");
my @mutant_type= <FILE>;
close(FILE);

open(FILE, "wild_BA");
my @wild_type= <FILE>;
close(FILE);

foreach (@mutant_type){

chomp ;
my @arr = split("\t", $_);

#Now wild type
foreach my $wild_line (@wild_type){
 chomp $wild_line;
my @wild_arr= split("\t", $wild_line);

if ($arr[0] eq $wild_arr[0]){
#print WRITE $a[1]."\n";

print 	WRITE "$wild_line\t$_\tSB\n";
	}

}
}



print "Wild_type_ID\tHLA\tWT_peptide\tWT_Score_EL\tWT_%Rank_EL\tWT_Score_BA\tWT_%Rank_BA\tWT_Aff(nM)\tWT_Bind_Level\tMutant_type_ID\tHLA\tMT_peptide\tMT_Score_EL\tMT_%Rank_EL\tMT_Score_BA\tMT_%Rank_BA\tMT_Aff(nM)\tMT_Bind_Level\tGene\tStart\tStop\tVariant_type\tHGVSc\tHGVSp\tProtein_Change\tVAF\tPolyphen_Effect\n";
close(WRITE);
#0_NEBL_208     NEBL    20809911        Missense        chr10                           10_NEBL_2081280 HLA-A*26:01     AAYKGVHPY       0.5522180       0.135   0.352645        0.376   10_NEBL_2081280 HLA-A*26:01     GVHPHIVEM       0.4734770       0.187   0.304858        0.584

#Find VAF
open(FILE, $ARGV[3]);
my @VAF = <FILE>;
close(FILE);

#READ THE SAMPLE_OUT
open(READ,"SAMPLE_OUT"); 

my @out =<READ>;
close(READ);

open(FILE, $ARGV[2]);
my @maf = <FILE>;
close(FILE);


my $FOUND=0;
my $vaf_found="";
my $vaf_effect="";
my $vaf_line="";
foreach my $line (@out){
chomp $line;


my @arr1 = split("\t", $line);
my ($gene, $pos, $NA) = split("_", $arr1[0]);

my @all = grep (/$gene/ , @VAF);


foreach (@maf){
chomp;

my @maf_arr= split("\t", $_);

if (($maf_arr[0] eq $gene) && ($maf_arr[5] =~ m/$pos/)){


foreach my $vaf (@all){
chomp $vaf;
my @vaf1 = split("\t", $vaf);

#print "$vaf1[0]\t$gene\n";

if (($vaf1[2] eq $gene) && ($vaf1[1] =~ m/$pos/)){

#print "$vaf1[2]\t$gene\t$vaf1[1]\t$pos\t1\n";
#FIND VAF and report otherwise 0 
#print "$line\t$maf_arr[0]\t$maf_arr[5]\t$maf_arr[6]\t$maf_arr[8]\t$maf_arr[34]\t$maf_arr[35]\t$maf_arr[54]\t$vaf1[5]\t$vaf1[6]\n";
$FOUND=1;
$vaf_found=$vaf1[5];
$vaf_effect =$vaf1[6];
$vaf_line = "$line\t$maf_arr[0]\t$maf_arr[5]\t$maf_arr[6]\t$maf_arr[8]\t$maf_arr[34]\t$maf_arr[35]\t$maf_arr[54]\t$vaf_found\t$vaf_effect";
}else{
#print "$vaf1[2]\t$gene\t$vaf1[1]\t$pos\t0\n";
#$FOUND =0;
#$vaf_found =0;
#$vaf_effect = "NA";
#$vaf_line= "$line\t$maf_arr[0]\t$maf_arr[5]\t$maf_arr[6]\t$maf_arr[8]\t$maf_arr[34]\t$maf_arr[35]\t$maf_arr[54]\t$vaf_found\t$vaf_effect";

#print "$vaf1[2]\t$gene\t$vaf1[1]\t$maf_arr[5]\n";
#print "$line\t$maf_arr[0]\t$maf_arr[5]\t$maf_arr[6]\t$maf_arr[8]\t$maf_arr[34]\t$maf_arr[35]\t$maf_arr[54]\t0\t$vaf1[6]\n";

}
}


if ($FOUND){

print "$vaf_line\n";
#print "$line\t$maf_arr[0]\t$maf_arr[5]\t$maf_arr[6]\t$maf_arr[8]\t$maf_arr[34]\t$maf_arr[35]\t$maf_arr[54]\t$vaf_found\t$vaf_effect\n";
}else{

print "$vaf_line\n";
#print "$line\t$maf_arr[0]\t$maf_arr[5]\t$maf_arr[6]\t$maf_arr[8]\t$maf_arr[34]\t$maf_arr[35]\t$maf_arr[54]\t$vaf_found\t$vaf_effect\n";
}



}
}
}








#print "Partial_id\tHugo_Symbol\tStart\tVariant_type\tChromosome\tNetMHCpan_ID\tWild_HLA-type\tWild_type_epitope\tWild_Score_EL\t%Wild_Rank_EL\tWild_Score_BA\t%Wild_Rank_BA\tMutant_NetMHCpan_ID\tMutant_HLA-type\t\tMutant_type_epitope\tMutant_Score_EL\t%Mutant_Rank_EL\tMutant_Score_BA\t%Mutant_Rank_BA\t";
