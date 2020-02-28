open F,"total.csv";
my $i;
while(<F>){
next if /sample/;
$i++;
my @line=split /,/;
next if $line[3] eq 1;
my $snp1=substr($line[4],0,1);
my $snp2=substr($line[4],1,1);
my $snp3=substr($line[5],0,1);
my $snp4=substr($line[5],1,1);
my $snp5=substr($line[6],0,1);
my $snp6=substr($line[6],1,1);
my $snp7=substr($line[7],0,1);
my $snp8=substr($line[7],1,1);
my $snp9=substr($line[8],0,1);
my $snp10=substr($line[8],1,1);
my @snp=($snp1,$snp2,$snp3,$snp4,$snp5,$snp6,$snp7,$snp8,$snp9,$snp10);
$line[3]=~s/0/1/;
$line[1]=~s/0/2/;
my $line2=join("\t",$i,$i,"0","0",$line[1],$line[3],@snp);
print "$line2\n";
}
