use strict;

my $Annotations= "BioMart_Ensembl100_Gencodev25_GRCm38p6_with_Transcript_Name.txt";
my $GeneList= "102522_Differential_Transcript_Expression_Significant_k=4.txt";
my $outfile="022723_Annotated_Differential_Transcript_Expression_Significant_k=4.txt";


open IN,"<$Annotations" or die $!;
my %Annotate;
while (<IN>){
  
  chomp;
  
  if (/ENSMUST(\d+)./){
    
    $Annotate{$1}=$_;
  }
}

close IN;

open OUT, ">$outfile" or die $!;
open IN2,"<$GeneList" or die $!;

while (<IN2>){
  
  if (/ENSMUST(\d+)./){
    
    my $ID=$1;
    
    my $match=$Annotate{$ID};
    if (defined $match){
      print "$ID\n";
      
      print OUT "$match$_";
      
    }
    
  }
}