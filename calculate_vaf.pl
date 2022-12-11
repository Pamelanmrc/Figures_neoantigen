use strict;
#perl calcualte_vaf.pl
open(FILE, $ARGV[0]);
while(<FILE>){
chomp;

if ($_ =~ m/^#/){

}else{


my @arr = split("\t", $_);
#print $arr[9]."\n";
my @missense = grep (/missense_variant/, $_);
#print @missense;
#print "$#missense\t$#arr\n";
if (scalar(@missense) <= 0){

}else{
my @mis = split("\t", $missense[0]);
#print @mis;
my @a = split(":", $mis[9]);
my @a1 = split(/\|/, $mis[7]);
#print "$arr[0]\t$arr[1]\t$a[2]\n";
print "$arr[0]\t$arr[1]\t$a1[4]\t$a1[7]\t$a1[5]\t$a[2]\t$a1[37]\n";

}
}

}
close(FILE);

