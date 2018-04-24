##This script will search for known TFBSs in the promoters of known genes from the Zea mays genome, count how many times each TFBS is found, and print this count to a file. 

#The next few lines are required in order to use Bio::Perl.
#!/usr/bin/perl
use lib "/opt/sharcnet/bioperl/1.6.923/lib/perl5";
use Bio::Perl;
use Bio::SeqIO;

#Use strict and use warnings forces Perl to provide more information on why a script wasn't executed. These warnings are very helpful and can save a lot of troubleshooting time.
use strict; 
use warnings;

#The script contains multiple print statements throughout that update the user as to what the script is currently working on. This is particularly useful since the script does take a fair amount of
#time to run.

#The next part of the code downloads and unzips the Zea mays GFF3 and genome FASTA files, which are required for the script to work. In order to do this, Perl uses the system() function which allows
#for the use of Unix commands within Perl. So, you would type the Unix command exactly how you would in the Unix terminal, and Perl will execute it as Unix would.

print("Now downloading and unzipping required files.\n");

system("wget ftp://ftp.ensemblgenomes.org/pub/plants/release-31/gff3/zea_mays/Zea_mays.AGPv3.31.gff3.gz");
system("gunzip Zea_mays.AGPv3.31.gff3.gz");
system("wget ftp://ftp.ensemblgenomes.org/pub/plants/release-31/fasta/zea_mays/dna/Zea_mays.AGPv3.31.dna.toplevel.fa.gz");
system("gunzip Zea_mays.AGPv3.31.dna.toplevel.fa.gz");

print("Done downloading and unzipping required files.\n");


#The next few lines use the system() function in order to use Unix commands as you would in the Unix terminal. The sed command is used to add ';' onto the end of every Maize gene name in the 
#maize.genes file so that when grep is used on the next code line, it will only retrieve exact matches for that gene name. The -i sed parameter is used to change the original input file (maize.genes) 
#so that an output file does not need to be specified. The -e sed parameter is used to tell sed that a regular expression is being used. When using grep, the -f parameter tells grep that it will be using
#patterns in a file (the gene names in the 'maize.genes') to extract lines from the maize GFF3 file. The results from this first grep are then piped into a second grep which then extracts lines that
#have 'ID=gene' in them. This is to make sure that, for each gene, we have the the earliest start position and the latest stop position (since each gene may have a record for multiple exons, which 
#each have a start and stop position within the gene itself, but doesn't give the complete picture). The results from this second grep are then written to a file called 'filteredGenes.txt'.
  
print("Now extracting the GFF3 line for each known Zea mays gene.\n");

system('sed -i -e "s/$/;/" maize.genes');
system("grep -f maize.genes Zea_mays.AGPv3.31.gff3 | grep 'ID=gene:' > filteredGenes.txt");

print("Done extracting the GFF3 line for each known Zea mays gene.\n");


print("Now extracting chromosome sequences for each chromosome in the Zea mays genome.\n");

#The next few lines use Bio::Perl to open the 'Zea_mays.AGPv3.31.dna.toplevel.fa' file as a FASTA file and save it to the '$seqIO' variable. This will later be used to create sequence objects for each
#record in the 'Zea_mays.AGPv3.31.dna.toplevel.fa' FASTA file.

my $seqdata="Zea_mays.AGPv3.31.dna.toplevel.fa"; #Will be changed to Zea_mays.AGPv3.31.dna.toplevel.fa
my $seqIO = new Bio::SeqIO(-file   => $seqdata,
                           -format => 'fasta');

my %chrom_seqs; #An empty hash that will later contain chromosome numbers as the keys and the sequence of the chromosome as the value. 

#The below while loop will iterate through each record in the $seqobject (which is the same as each record in the 'Zea_mays.AGPv3.31.dna.toplevel.fa' file), and split the $seqobj at each space, saving 
#the record ID to the '@name' array, and saving the first field of the '@name' element, which is the chromosome number, to the '$chrom' variable. The loop will then define a key:value pair
#for the current $seqobj, where the key is the chromosome number and the value is the chromosome sequence and save this key:value pair in the '%chrom_seqs' hash.

while( my $seqobj = $seqIO->next_seq ) {
    my @name = split / /, $seqobj->display_id();
    my $chrom = shift @name;
    $chrom_seqs{$chrom} = $seqobj->seq();
}

print("Done extracting chromosome sequences for each chromosome in the Zea mays genome.\n");


print("Now extracting promoter sequences for each Zea mays gene.\n");

my @genePromoters = (); #An empty array that will eventually hold the promoter sequence for each gene in the '@filteredGenes' array.

#The following foreach loop will go through each line in the 'filteredGenes' file one at a time. The split() function will then split the current line at each tab character (\t) and save each 
#piece to the '@filteredGeneSplit' array. The loop then finds the chromosome number (first element in the '@filteredGeneSplit' array), gene start position (third element in '@filteredGeneSplit'),
#the gene stop position (fourth element in '@filteredGeneSplit'), and the strand orientation (sixth position in '@filteredGeneSplit'). The chromosome number is required to find the relevant chromosome
#sequence for the current gene. The start position is important to find the promoter for genes encoded on the '+' strand ($promoter = substr $chrom_seqs{$gff3chrom}, $startPosition-500-1, 500;).
#The stop position is important for finding the promoter for genes encoded on the '-' strand ($promoter = substr $chrom_seqs{$gff3chrom}, $stopPosition+500+1, 500; 
#$revcomp_Promoter = revcom_as_string($promoter);. The revcom_as_string() function finds the reverse complement of $promoter, which is required for genes encoded on the '-' strand). Strand orientation
#is important for knowing which gene is encoded on which strand. An if/else statement is then used to tell Perl how to find the promoters for genes on the '+' or '-' strand, since as written previously,
#the approaches are different. The ($promoter =~ s/^.*N//i) line makes Perl keep the sequence after the ^ and up to the first occurrence of an N (or n, since \i makes the match case insensitive) 
#so that any low quality part of the sequence is not used. 

open (my $filteredGenes, '<', 'filteredGenes.txt') or die "Could not open 'filteredGenes.txt' for reading: $!\n";
while (my $filteredGeneLine = <$filteredGenes>) {
	my @filteredGeneSplit = split /\t/, $filteredGeneLine;
        my $gff3chrom = $filteredGeneSplit[0];
        my $startPosition = $filteredGeneSplit[3];
	my $stopPosition = $filteredGeneSplit[4];
        my $strandOrientation = $filteredGeneSplit[6];
	if ($strandOrientation eq '+') {
		my $promoter = substr $chrom_seqs{$gff3chrom}, $startPosition-500-1, 500;
		$promoter =~ s/^.*N//i; 
		push @genePromoters, $promoter;
	} else {
		my $promoter = substr $chrom_seqs{$gff3chrom}, $stopPosition, 500;
		my $revcomp_Promoter = revcom_as_string($promoter);
		$revcomp_Promoter =~ s/^.*N//i;
		push @genePromoters, $revcomp_Promoter;
	}
}

print("Done extracting promoter sequences for each Maize gene.\n");


print("Now finding TFBSs present in each promoter.\n");

my @matchedTFBS = (); #An empty array that will eventually hold the sequence of any matched TFBS for each promoter in the '@genePromoters' array.

#The following foreach loop iterates through each promoter in the '@genePromoters' array one at a time. The foreach loop will also open the 'promoters' file for reading (or exit the program if this file
#cannot be opened for some reason and print a message telling the user why the program has been ended prematurely) and close the 'promoters' file with each iteration so that Perl will start at the 
#beginning of this file for each promoter in the '@genePromoters' array. A while loop then iterates through each line (or TFBS) in the 'promoters' file, and chomps each line so that any newline
#characters do not alter any results. An if statement is then used to determine if the current TFBS is found in the current promoter sequence (if ($maizePromoters =~ m/$maizeTFBSLines/i) - the
#/i designation is used to make the regex is case insensitive, since some TFBS use lower-case letters but the promoter sequences are typically upper-case. If a match is found, then the TFBS, or
#$maizeTFBSLines line, will be pushed onto the '@matchedTFBS' array.

foreach my $maizePromoters (@genePromoters) {
        open (my $maizeTFBS, '<', 'promoters') or die "Could not open 'promoters' for reading: $!\n";
        while (my $maizeTFBSLines = <$maizeTFBS>) {
                chomp($maizeTFBSLines);
		if ($maizePromoters =~ m/$maizeTFBSLines/i) {
			push @matchedTFBS, $maizeTFBSLines;
		}
	}
	close($maizeTFBS);
}

print("Done finding TFBSs present in each promoter.\n");


print("Now determining the number of times each TFBS has been found in the maize promoters.\n");

#The following lines open the 'promoters' file for reading and the 'Count_Of_Maize_TFBSs_Found.txt' file for writting. If for some reason these files cannot be opened, then the or die() function
#will stop the program, and print a message to the user as to why the program has been stopped. 

open (my $maizeTFBSs, '<', 'promoters') or die "Could not open 'promoters' for reading: $!\n";
open (my $TFBScount, '>', 'Count_Of_Maize_TFBSs_Found.txt') or die "Could not open 'Count_Of_Maize_TFBS_Found.txt' for writting: $!\n";

#The following while loop iterates through each TFBS in the 'promoter' file, and chomps each line so that newline characters do not affect any results. Within the while loop, a 
#'$counter' variable is set to 0. This counter will keep track of how many times each TFBS occurs in the '@matchedTFBS' array. A foreach loop will go through each element in the '@matchedTFBS' 
#array, and if the current element matches the current TFBS, then the counter will be increased by one. After every element in the '@matchedTFBS' array has been checked, the TFBS sequence 
#and number of times the TFBS has occured will be written to the 'Count_Of_Maize_TFBSs_Found.txt' file. 

while (my $maizeTFBSLines = <$maizeTFBSs>) {
	chomp($maizeTFBSLines);
	my $counter = 0;
	foreach my $matchedTFBSs (@matchedTFBS) {
		if ($matchedTFBSs eq $maizeTFBSLines) {
			$counter ++;
		}
	}
	print $TFBScount ("The TFBS $maizeTFBSLines occurs $counter time(s).\n");
}

print("Done counting the number of times each TFBS has been found in the Maize promoters.\n");
print("A report containing each TFBS sequence and the number of times it was found in the Zea mays promoters can be found in the 'Count_Of_Maize_TFBSs_Found.txt' file.\n");
print("You can find this file in the same directory that you ran the Perl script from.\n");


