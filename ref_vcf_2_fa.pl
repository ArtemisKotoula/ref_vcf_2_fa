#!/usr/bin/perl

#######################################################
# The script reconstructs the coding sequence of each sample
# by introducing the SNPs from the VCF into the reference sequence.

#######################################################
# PATHS for vcf and reference files
$file="vcf.txt"; #Path to vcf file
$fasta="one_line.fa";   #Path to a "oneline" FASTA file that contains the
                        # reference sequence
############################################################

use File::Path qw(make_path);

# Reads the vcf file 
open (FIN, "< $file") || die "Could not open vcf file\n$!\n";
#make a list with each trans id
while (my $line=<FIN>) {
    chomp($line);
    #ignores lines that start with #
    if ($line !~ /^#/){
        #print STDOUT "$line\n";
        @line=split(/\t/, $line);   #splits line with tabs
        $trans_id=$line[3];#
        push(@list_of_trans_ids,$trans_id);
        
        $gene_name=$line[6];
        
        push(@list_of_gene_names, $gene_name);
        #get the samples so 
        #$samples=$line[5];
    }
}
close(FIN);

#print STDOUT "@list_of_gene_names\n";


$samples=$line[5];
my $char = "/";
my $number_of_samples = () = $samples =~ /\Q$char/g;
#print "count<$count> of <$char> in <$samples>\n";
print STDOUT "Number of samples detected: $number_of_samples\n";

#keep uniq trans ids
   my %seen = ();
   my @unique_trans_ids = grep { ! $seen{ $_ }++ } @list_of_trans_ids;
   @unique_gene_names = grep { ! $seen{ $_ }++ } @list_of_gene_names;


#print STDOUT @unique_gene_names;

$samples="";
$length=1+$number_of_samples*2;

# For each transcript:
# - searches the FASTA
# - finds the reference CDS
# - saves it as $ref_seq

make_path("./ref_out") unless -d "ref_out";

$name_index=0;
foreach $elem (@unique_trans_ids){
    @lines_per_sample=();
    @seq_per_sample=();
    @each_sample=();

    open(FASTA, "< $fasta") || die "Could not open fasta file\n";
        #loops to find transcript id
        while ($fasta_header=<FASTA>){
            chomp($fasta_header);
            $fasta_seq=<FASTA>; #REFERENCE sequence
            chomp($fasta_seq);
            if ($fasta_header =~ $elem){
                $ref_seq=$fasta_seq;
            }
        }
    close(FASTA);

        #reference sequenves to create, for each sample
        for ($i=0; $i < $length; $i++ ){
            push(@seq_per_sample,$ref_seq);
        }

        #print STDOUT "seq per sample length:$length\n";

    open (FIN, "< $file") || die "Could not open vcf file\n$!\n";
    while (my $line=<FIN>) {
        chomp($line);
        if (($line =~ $elem) && ($line !~ /^#/)){

            @line=split(/\t/, $line);   #splits line with tabs
            $trans_id=$line[3];#
            $replace=$line[4];
            $samples=$line[5];
            $gene_name=$line[6];
            $gene_id=$line[7];
            $positions=$line[8];

            #find al1:
            my @replace=split(/>/, $replace);
            $al1=$replace[-1];#
            
            #print STDOUT "$al1\n";

            #find position:
            my @positions=split(/\//, $positions);
            $CDS_pos=$positions[0];#
            $array_CDS_POS=$CDS_pos-1;

            #samples
            $samples =~ s/\///g;
            @each_sample=split(//, $samples);

            for ($i=0; $i < $length; $i++){
                $sequence=$seq_per_sample[$i];
                @seq=split(//, $sequence); #split sequence

                if ($each_sample[$i] eq "1"){
                    $seq[$array_CDS_POS] = $al1;

                } elsif ($each_sample[$i] eq "0"){
                    #no change
                } elsif ($each_sample[$i] eq ".") {
                    $seq[$array_CDS_POS] = "N";
                }

                $sequence=join("", @seq);
                $seq_per_sample[$i]=$sequence;

            }

        }    

    }
    close(FIN);

    open(OUT, "> ./ref_out/$unique_gene_names[$name_index].$elem.fa") or die "Cannot write output file: $!";
    
    $index=1;
    for ($i=0; $i < $length; $i=$i+2){
        $number=$i+1;
        
        if ($i == ($length -1)) {
            print OUT ">$elem"."_reference_seq\n$seq_per_sample[$i]\n";
        } else{
            $sample_name="sample".$index;
        print OUT ">$sample_name\ts1\t$unique_gene_names[$name_index]\t$elem\n$seq_per_sample[$i]\n";
        print OUT ">$sample_name\ts2\t$unique_gene_names[$name_index]\t$elem\n$seq_per_sample[$number]\n";
        }
        $index++;
    }
    close(OUT);
    $name_index++;
}
