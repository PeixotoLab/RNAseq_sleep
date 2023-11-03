use strict;

my $Annotations= "BioMart_Ensembl100_Gencodev25_GRCm38p6.txt";
my $GeneList= "DGE_Significant_k=4.txt";
my $outfile="DGE_Significant_k=4_Annotated.txt";


open IN,"<$Annotations" or die $!;
my %Annotate;
while (<IN>){
  
  chomp;
  
  if (/ENSMUSG(\d+)./){
    
    $Annotate{$1}=$_;
  }
}

close IN;

open OUT, ">$outfile" or die $!;
open IN2,"<$GeneList" or die $!;

while (<IN2>){
  
  if (/ENSMUSG(\d+)./){
    
    my $ID=$1;
    
    my $match=$Annotate{$ID};
    if (defined $match){
      print "$ID\n";
      
      print OUT "$match$_";
      
    }
    
  }
}
