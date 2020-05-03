#!/usr/bin/env perl

=head
Elloreas, a genome assembler. https://github.com/shelkmike/Elloreas
=cut
$Elloreas_version="1.15";


#determining the path to the folder where Elloreas is located. It is needed to run auxiliary scripts like plot_coverage_for_Elloreas.r
use Cwd 'abs_path';
$path_to_the_folder_where_Elloreas_is_located=abs_path($0);
$path_to_the_folder_where_Elloreas_is_located=~s/^(.+\/)\S.+/$1/; #removing the script name, retaining only the path to the folder

&take_parameters_from_the_command_line; #the subroutine to take parameters provided by the user in a command line

#Creating the output folder. Elloreas doesn't remove the previous folder if it exists already, because otherwise if a user provides a wrong output folder name, Elloreas may delete important data.
system("mkdir --parents $output_folder");
open LOGFILE, "> $output_folder/elloreas_logs.txt";


&check_whether_required_parameters_and_scripts_exist; #this subroutine checks whether programs required by Elloreas exist in $PATH and whether there are scripts plot_coverage_for_Elloreas.r and make_dotplot_for_Elloreas.r in the same folder where elloreas.pl lies.
&print_parameters_provided_by_the_user_via_the_command_line; #this subroutine prints a list of parameters provided by the user via the command line.
&check_that_the_input_FASTA_file_contains_only_one_sequence; #checking that the input file contains only one sequence.



open HISTORY_OF_ALTERNATIVE_EXTENSIONS, "> $output_folder/history_of_alternative_extensions.txt"; #the file where I write top-5 alternative extensions at each iteration (or less, if there are less than five).

#loading the starter sequence
open CONTIG_FILE, "< $starter_file_path";
$contig_sequence="";
$contig_title="";
while(<CONTIG_FILE>)
{
	if($_=~/^>(.+)$/)
	{
		$contig_title=$1;
	}
	if($_!~/^>/)
	{
		$contig_sequence.=$_;
		chomp($contig_sequence);
	}
}
#removing space characters in the sequence of the starter
$contig_sequence=~s/\s//g;
#converting the starter sequence to the upper case
$contig_sequence=uc($contig_sequence);
close(CONTIG_FILE);

#chopping $number_of_bases_to_trim_from_the_starter_edges_in_the_very_beginning_of_the_assembly bases from each end. If the sequence is so short that after this chopping it becomes shorter than , then die.
if(length($contig_sequence)>=(2*$number_of_bases_to_trim_from_the_starter_edges_in_the_very_beginning_of_the_assembly+$min_aligned_length))
{
	$contig_sequence=substr($contig_sequence,$number_of_bases_to_trim_from_the_starter_edges_in_the_very_beginning_of_the_assembly,length($contig_sequence)-2*$number_of_bases_to_trim_from_the_starter_edges_in_the_very_beginning_of_the_assembly);
	print LOGILE "The contig length was ".($contig_sequence+2*$number_of_bases_to_trim_from_the_starter_edges_in_the_very_beginning_of_the_assembly).". After chopping $number_of_bases_to_trim_from_the_starter_edges_in_the_very_beginning_of_the_assembly bases from each edge according to the value of the --number_of_bases_to_trim_from_the_starter_edges parameter it became ".length($contig_sequence)." bases long.";
}
else
{
	print LOGFILE "\nEnding, because the starter is too short. You may want to reduce the value of the \"--number_of_bases_to_trim_from_the_starter_ends\" parameter\n";
	print HISTORY_OF_ALTERNATIVE_EXTENSIONS "\nElloreas has finished. See elloreas_logs.txt\n";
	open CONTIG_FROM_THE_FINAL_ITERATION, "> $output_folder/contig_from_the_final_iteration.fasta";
	print CONTIG_FROM_THE_FINAL_ITERATION ">$contig_title\n";
	print CONTIG_FROM_THE_FINAL_ITERATION "$contig_sequence\n";
	close(CONTIG_FROM_THE_FINAL_ITERATION);
	&make_the_final_coverage_plot_and_dot_plot();
	exit();
}

#now, starting the iterations
foreach $iteration_number(1..$max_number_of_cycles)
{
	$current_edge_sequence="";
	#if the user asked for this, remove data from all previous iterations
	if($should_elloreas_delete_folders_of_all_iterations=~/true/)
	{
		system("rm -rf ./$output_folder/Iteration*");
	}
	system("mkdir $output_folder/Iteration$iteration_number");
	#taking the right edge of the contig to map reads
	#if the contig is shorter than $contig_edge_size_to_use_for_alignment then we take the whole contig. If it is longer, take only the right $contig_edge_size_to_use_for_alignment bases
	if(length($contig_sequence)>=$contig_edge_size_to_use_for_alignment)
	{
		$current_edge_sequence=substr($contig_sequence,length($contig_sequence)-$contig_edge_size_to_use_for_alignment, $contig_edge_size_to_use_for_alignment);
		#creating the file with the edge's sequence
		open CURRENT_EDGE_OUTFILE, "> $output_folder/Iteration$iteration_number/Edge_iteration$iteration_number.fasta";
		print CURRENT_EDGE_OUTFILE ">Edge_iteration$iteration_number\n$current_edge_sequence\n";
		close(CURRENT_EDGE_OUTFILE);
	}
	else
	{
		$current_edge_sequence=$contig_sequence;
		#creating the file with the edge's sequence
		open CURRENT_EDGE_OUTFILE, "> $output_folder/Iteration$iteration_number/Edge_iteration$iteration_number.fasta";
		print CURRENT_EDGE_OUTFILE ">Edge_iteration$iteration_number\n$current_edge_sequence\n";
		close(CURRENT_EDGE_OUTFILE);
	}
	print LOGFILE "=============================================================\n";
	print LOGFILE "Iteration $iteration_number\n";
	print LOGFILE "=============================================================\n";
	print HISTORY_OF_ALTERNATIVE_EXTENSIONS "=============================================================\n";
	print HISTORY_OF_ALTERNATIVE_EXTENSIONS "Iteration $iteration_number\n";
	print HISTORY_OF_ALTERNATIVE_EXTENSIONS "=============================================================\n";

	
	@array_extensions=(); #array with parts of reads overhanging the edge
	$average_length_of_extensions=0; #average length of parts of reads overhanging the edge
	
	#mapping reads
	system("minimap2 -x $sequencing_technology -a --secondary no --end-bonus 100 --sam-hit-only -o $output_folder/Iteration$iteration_number/mapping.it$iteration_number.sam -t $cores_to_use $output_folder/Iteration$iteration_number/Edge_iteration$iteration_number.fasta $reads_file_path");
	
	#if required, make a coverage plot. To do this, I first need to convert sam to bam
	if($should_elloreas_draw_coverage_plots=~/true/)
	{
		system("samtools view -Sbh $output_folder/Iteration$iteration_number/mapping.it$iteration_number.sam >$output_folder/Iteration$iteration_number/mapping.it$iteration_number.bam");
		system("samtools sort $output_folder/Iteration$iteration_number/mapping.it$iteration_number.bam >$output_folder/Iteration$iteration_number/mapping.it$iteration_number.sorted.bam");
		system("samtools index $output_folder/Iteration$iteration_number/mapping.it$iteration_number.sorted.bam");
		system("bedtools genomecov -d -ibam $output_folder/Iteration$iteration_number/mapping.it$iteration_number.sorted.bam > $output_folder/Iteration$iteration_number/coverage_list.it$iteration_number.txt");
		system("Rscript $path_to_the_folder_where_Elloreas_is_located/plot_coverage_for_Elloreas.r $output_folder/Iteration$iteration_number/coverage_list.it$iteration_number.txt $output_folder/Iteration$iteration_number/coverage_plot.it$iteration_number.jpeg");
	}
	
	#first, I go through the sam-file to determine the median length of read's overhangs. (but counting it only using reads for which at least $min_aligned_length bases were mapped)
	@array_of_overhang_lengths=();
	$median_overhang_length=0;
	open SAMFILE, "< $output_folder/Iteration$iteration_number/mapping.it$iteration_number.sam";               
	while(<SAMFILE>)
	{
	#If I see a string like '@SQ     SN:contig_1-2-3 LN:152189', I take the reference sequence length from it. It is needed to understand which reads overlap the edge.
		if($_=~/^\@SQ\s+SN:(.*?)\s+LN:(\d+)/)
		{
			$sequence_we_mapped_to_title=$1;
			$sequence_we_mapped_to_length=$2;
		}
		#для каждого рида я буду мерять его длину с помощью length(последовательность). Если координата начала рида (а это его левая координата, по стандарту формата sam) плюс его длина больше, чем длина последовательности, на которую мы картировали, то этот рид нам подходит. I also take into account deletions in the read (\d+D) and insertions in the read (\d+I).

		
		#m64035_191023_221642/118/ccs    0       Capsella_bursa_pastoris_Msk_mitochondrion       214903  60      2241S359M1D150M4480S    *       0       0       CGAAGAAAGGCATGCGAGAAAAGCATATTGGCTAGTGATTGTGAGGCTCCAATTCTTGACTGGAGTGGACACCAAAGGCCTCCGCCCCCCCATCCTTTGGATAGAAGGGCAGAATTTTAGGTTTTTTCATGTTGTCAAAGAGTTGAACAATGGTTTTTTCGTGTTGTCAAAGAGTTGAACTATGAAAATAGATGGCGAGTGCCTGATCGAATTGATCAGGTCATGTAGGAACAAGGTTCAAGTCTACCGGTCTGTTAGGATGCCTCAGCTGC...
		
		#если это строка с картированием рида
		if($_=~/^\S+\s+\S+\s+$sequence_we_mapped_to_title\s+/)
		{
			@string_split=split(/\s+/, $_);
			$start_position_of_the_mapped_part_of_the_read=$string_split[3]; #this is the start position of the MAPPED part of the read, i.e. not counting left hard- and soft-clipped regions
			$read_sequence=$string_split[9];
			$CIGAR_string=$string_split[5];
			
			#calculate the number of differences between a read and the reference in the mapped part of the read.
			if($_=~/\tNM:i:(\d+)\t/)
			{
				$number_of_differences_in_the_mapped_part=$1;
			}
			else
			{
				print LOGFILE "died because cannot find the number of differences (\"NM:i:\\d+\") in the SAM file";
				die "died because cannot find the number of differences (\"NM:i:\\d+\") in the SAM file";
			}
			
			$temp_CIGAR_string=$CIGAR_string;
			#I remove the information about hard clipping from CIGAR, because it isn't important for anything.
			$temp_CIGAR_string=~s/\d+H//g;
			#If there was soft clipping at the left side of the read, I remove that number of bases from the read's left end. Because the start position of the read is indicated only for the mapped part, while the sequence that is soft-clipped is included in the read sequence in the 10th column ($string_split[9]) of the SAM file
			if($temp_CIGAR_string=~/^(\d+)S/)
			{
				$number_of_soft_clipped_bases_at_the_left_side_of_the_read=$1;
				$temp_CIGAR_string=~s/^\d+S//;
				$read_sequence=substr($read_sequence,$number_of_soft_clipped_bases_at_the_left_side_of_the_read,$read_sequence-$number_of_soft_clipped_bases_at_the_left_side_of_the_read);
			}
			
			#counting the number of deletions in the read (\d+D) and insertions in the read (\d+I).

			$number_of_I_letters_in_CIGAR=0;
			$number_of_D_letters_in_CIGAR=0;
			$number_of_S_letters_at_the_right_side_of_CIGAR=0;
			$number_of_M_letters_in_CIGAR=0;
			
			while($temp_CIGAR_string=~/(\d+)I/)
			{
				$number_of_I_letters_in_CIGAR+=$1;
				$temp_CIGAR_string=~s/\d+I//;
			}
			while($temp_CIGAR_string=~/(\d+)D/)
			{
				$number_of_D_letters_in_CIGAR+=$1;
				$temp_CIGAR_string=~s/\d+D//;
			}
			if($temp_CIGAR_string=~/(\d+)S$/)
			{
				$number_of_S_letters_at_the_right_side_of_CIGAR=$1;
			}
			while($temp_CIGAR_string=~/(\d+)M/)
			{
				$number_of_M_letters_in_CIGAR+=$1;
				$temp_CIGAR_string=~s/\d+M//;
			}
			
						
			#calculating the sequence similarity between the mapped part of the read and the contig. Indels are already counted in NM:i:(\d+)
			$sequence_similarity_between_the_mapped_part_of_the_read_and_the_contig=($number_of_M_letters_in_CIGAR-$number_of_differences_in_the_mapped_part)/$number_of_M_letters_in_CIGAR;
			
			$read_length=length($read_sequence);
			
			$number_of_deleted_bases_in_the_read_minus_the_number_of_inserted_bases_in_the_read=$number_of_D_letters_in_CIGAR-$number_of_I_letters_in_CIGAR;
			
			#checking that the right does not have an unmapped part before the edge of the contig (i.e. the hanging part starts exactly at the contig edge and not earlier. It can start earlier if the region of a read that is adjacent (from the left side) to the contig's edge contains so many sequencing errors that it does no align). This is needed because of the problem described in the description for the version 1.10, point "2)".
			$does_the_read_has_an_unmapped_part_before_the_edge_of_the_contig="no"; #may be "no" or "yes"
			if(($read_length-$number_of_S_letters_at_the_right_side_of_CIGAR+$number_of_deleted_bases_in_the_read_minus_the_number_of_inserted_bases_in_the_read)!=(length($current_edge_sequence)-$start_position_of_the_mapped_part_of_the_read+1))
			{
				$does_the_read_has_an_unmapped_part_before_the_edge_of_the_contig="yes";
				#print LOGFILE "this read has its part, which is close to the contig edge, unmapped\n"; 
			}
			else
			{
				#print LOGFILE "this read has its part, which is close to the contig edge, mapped\n";
			}
			
			if((($start_position_of_the_mapped_part_of_the_read+$read_length+$number_of_deleted_bases_in_the_read_minus_the_number_of_inserted_bases_in_the_read)>$sequence_we_mapped_to_length)&&($number_of_M_letters_in_CIGAR>=$min_aligned_length)&&($sequence_similarity_between_the_mapped_part_of_the_read_and_the_contig>=$min_sequence_similarity)&&($does_the_read_has_an_unmapped_part_before_the_edge_of_the_contig!~/yes/)&&($number_of_S_letters_at_the_right_side_of_CIGAR>0))
			{
				$overhang_length_of_this_read=$number_of_S_letters_at_the_right_side_of_CIGAR;
				push(@array_of_overhang_lengths,$overhang_length_of_this_read);
			}
		}

	}
	close(SAMFILE);
	
	#calculating the median overhang length. If the number of elements is even, I take the smaller of the two middle values.
	@sorted_array_of_overhang_lengths=sort {$a <=> $b} @array_of_overhang_lengths;
    $number_of_elements=$#sorted_array_of_overhang_lengths+1;
	$rank_of_the_middle_element=int($number_of_elements/2); #rank is being counted from zero here, not from 1.
	$median_overhang_length=$sorted_array_of_overhang_lengths[$rank_of_the_middle_element];
	
	print LOGFILE "The length of the contig to be extended is ".length($contig_sequence)." bases\n";
	if($median_overhang_length!~/^$/) #if there is at least one overhang
	{
		print LOGFILE "The median overhang length is $median_overhang_length bases. There are ".($rank_of_the_middle_element+1)." reads with overhangs this long or longer.\n"; #($rank_of_the_middle_element+1) because the rank was counted from zero, not from 1.
	}
	else
	{
		#if there are no overhangs Elloreas will write nothing, because anyway it will finish soon and write "Ending, because there are no extensions".
	}
	
	#теперь открываю файл с картированием и нахожу там все риды, картировавшиеся на край. Strolling through the same sam-file twice doesn't take much time, because I have directed minimap2 to print only mapped reads.
	
	open SAMFILE, "< $output_folder/Iteration$iteration_number/mapping.it$iteration_number.sam";
	$number_of_overhanging_reads_that_map_well_enough=0; #the number of reads that map to the edge with enough sequence similarity and long enough mapping part.
	while(<SAMFILE>)
	{
	#If I see a string like '@SQ     SN:contig_1-2-3 LN:152189', I take the reference sequence length from it. It is needed to understand which reads overlap the edge.
		if($_=~/^\@SQ\s+SN:(.*?)\s+LN:(\d+)/)
		{
			$sequence_we_mapped_to_title=$1;
			$sequence_we_mapped_to_length=$2;
		}
	
		#для каждого рида я буду мерять его длину с помощью length(последовательность). Если координата начала рида (а это его левая координата, по стандарту формата sam) плюс его длина больше, чем длина последовательности, на которую мы картировали, то этот рид нам подходит. I also take into account deletions in the reference (\d+D) and insertions in the reference (\d+I). Of course, in fact they are, respectively, insertions and deletions in reads, by the SAM specification calls them so (https://samtools.github.io/hts-specs/SAMv1.pdf).

		
		#m64035_191023_221642/118/ccs    0       Capsella_bursa_pastoris_Msk_mitochondrion       214903  60      2241S359M1D150M4480S    *       0       0       CGAAGAAAGGCATGCGAGAAAAGCATATTGGCTAGTGATTGTGAGGCTCCAATTCTTGACTGGAGTGGACACCAAAGGCCTCCGCCCCCCCATCCTTTGGATAGAAGGGCAGAATTTTAGGTTTTTTCATGTTGTCAAAGAGTTGAACAATGGTTTTTTCGTGTTGTCAAAGAGTTGAACTATGAAAATAGATGGCGAGTGCCTGATCGAATTGATCAGGTCATGTAGGAACAAGGTTCAAGTCTACCGGTCTGTTAGGATGCCTCAGCTGC...
		
		#если это строка с картированием рида
		if($_=~/^\S+\s+\S+\s+$sequence_we_mapped_to_title\s+/)
		{
			@string_split=split(/\s+/, $_);
			$start_position_of_the_mapped_part_of_the_read=$string_split[3]; #this is the start position of the MAPPED part of the read, i.e. not counting left hard- and soft-clipped regions
			$read_sequence=$string_split[9];
			$CIGAR_string=$string_split[5];
			
			#calculate the number of differences between a read and the reference in the mapped part of the read.
			if($_=~/\tNM:i:(\d+)\t/)
			{
				$number_of_differences_in_the_mapped_part=$1;
			}
			else
			{
				print LOGFILE "died because cannot find the number of differences (\"NM:i:\\d+\") in the SAM file";
				die "died because cannot find the number of differences (\"NM:i:\\d+\") in the SAM file";
			}
			
			$temp_CIGAR_string=$CIGAR_string;
			#I remove the information about hard clipping from CIGAR, because it isn't important for anything.
			$temp_CIGAR_string=~s/\d+H//g;
			#If there was soft clipping at the left side of the read, I remove that number of bases from the read's left end. Because the start position of the read is indicated only for the mapped part, while the sequence that is soft-clipped is included in the read sequence in the 10th column ($string_split[9]) of the SAM file
			if($temp_CIGAR_string=~/^(\d+)S/)
			{
				$number_of_soft_clipped_bases_at_the_left_side_of_the_read=$1;
				$temp_CIGAR_string=~s/^\d+S//;
				$read_sequence=substr($read_sequence,$number_of_soft_clipped_bases_at_the_left_side_of_the_read,$read_sequence-$number_of_soft_clipped_bases_at_the_left_side_of_the_read);
			}
			
			#counting the number of deletions in the read (\d+D) and insertions in the read (\d+I).

			$number_of_I_letters_in_CIGAR=0;
			$number_of_D_letters_in_CIGAR=0;
			$number_of_S_letters_at_the_right_side_of_CIGAR=0;
			$number_of_M_letters_in_CIGAR=0;
			
			while($temp_CIGAR_string=~/(\d+)I/)
			{
				$number_of_I_letters_in_CIGAR+=$1;
				$temp_CIGAR_string=~s/\d+I//;
			}
			while($temp_CIGAR_string=~/(\d+)D/)
			{
				$number_of_D_letters_in_CIGAR+=$1;
				$temp_CIGAR_string=~s/\d+D//;
			}
			if($temp_CIGAR_string=~/(\d+)S$/)
			{
				$number_of_S_letters_at_the_right_side_of_CIGAR=$1;
			}
			while($temp_CIGAR_string=~/(\d+)M/)
			{
				$number_of_M_letters_in_CIGAR+=$1;
				$temp_CIGAR_string=~s/\d+M//;
			}
			
			#calculating the sequence similarity between the mapped part of the read and the contig. Indels are already counted in NM:i:(\d+)
			$sequence_similarity_between_the_mapped_part_of_the_read_and_the_contig=($number_of_M_letters_in_CIGAR-$number_of_differences_in_the_mapped_part)/$number_of_M_letters_in_CIGAR;
			
			#print "\$sequence_similarity_between_the_mapped_part_of_the_read_and_the_contig is $sequence_similarity_between_the_mapped_part_of_the_read_and_the_contig\n";
			
			$read_length=length($read_sequence);
			
			$number_of_deleted_bases_in_the_read_minus_the_number_of_inserted_bases_in_the_read=$number_of_D_letters_in_CIGAR-$number_of_I_letters_in_CIGAR;
			
			#print "CIGAR is $CIGAR_string . Reduced CIGAR is $temp_CIGAR_string . \$number_of_deleted_bases_in_the_read_minus_the_number_of_inserted_bases_in_the_read is $number_of_deleted_bases_in_the_read_minus_the_number_of_inserted_bases_in_the_read\n";
			
			#checking that the right does not have an unmapped part before the edge of the contig (i.e. the hanging part starts exactly at the contig edge and not earlier. It can start earlier if the region of a read that is adjacent (from the left side) to the contig's edge contains so many sequencing errors that it does no align). This is needed because of the problem described in the description for the version 1.10, point "2)".
			$does_the_read_has_an_unmapped_part_before_the_edge_of_the_contig="no"; #may be "no" or "yes"
			if(($read_length-$number_of_S_letters_at_the_right_side_of_CIGAR+$number_of_deleted_bases_in_the_read_minus_the_number_of_inserted_bases_in_the_read)!=(length($current_edge_sequence)-$start_position_of_the_mapped_part_of_the_read+1))
			{
				$does_the_read_has_an_unmapped_part_before_the_edge_of_the_contig="yes";
				#print LOGFILE "this read has its part, which is close to the contig edge, unmapped\n"; 
			}
			else
			{
				#print LOGFILE "this read has its part, which is close to the contig edge, mapped\n";
			}
			
			#printing to LOGS, if the part close to the edge is not mapped (this should be commented in the future)
			if((($start_position_of_the_mapped_part_of_the_read+$read_length+$number_of_deleted_bases_in_the_read_minus_the_number_of_inserted_bases_in_the_read)>$sequence_we_mapped_to_length)&&($number_of_M_letters_in_CIGAR>=$min_aligned_length)&&($sequence_similarity_between_the_mapped_part_of_the_read_and_the_contig>=$min_sequence_similarity)&&($does_the_read_has_an_unmapped_part_before_the_edge_of_the_contig!~/yes/))
			{
				#print LOGFILE "this read has its part, which is close to the contig edge, mapped. The CIGAR string is $CIGAR_string . \n";
			}
			elsif((($start_position_of_the_mapped_part_of_the_read+$read_length-$number_of_deleted_bases_in_the_read_minus_the_number_of_inserted_bases_in_the_read)>$sequence_we_mapped_to_length)&&($number_of_M_letters_in_CIGAR>=$min_aligned_length)&&($sequence_similarity_between_the_mapped_part_of_the_read_and_the_contig>=$min_sequence_similarity)&&($does_the_read_has_an_unmapped_part_before_the_edge_of_the_contig=~/yes/))
			{
				#print LOGFILE "this read has its part, which is close to the contig edge, unmapped. ($read_length-$number_of_S_letters_at_the_right_side_of_CIGAR+$number_of_deleted_bases_in_the_read_minus_the_number_of_inserted_bases_in_the_read)!=(".length($current_edge_sequence)."-$start_position_of_the_mapped_part_of_the_read+1) . The CIGAR string is $CIGAR_string . \n"; 
			}
			
			
			if((($start_position_of_the_mapped_part_of_the_read+$read_length-$number_of_deleted_bases_in_the_read_minus_the_number_of_inserted_bases_in_the_read)>$sequence_we_mapped_to_length)&&($number_of_M_letters_in_CIGAR>=$min_aligned_length)&&($sequence_similarity_between_the_mapped_part_of_the_read_and_the_contig>=$min_sequence_similarity)&&($does_the_read_has_an_unmapped_part_before_the_edge_of_the_contig!~/yes/)&&($number_of_S_letters_at_the_right_side_of_CIGAR>0))
			{

				$number_of_overhanging_reads_that_map_well_enough++;
				push(@array_extensions, substr($read_sequence, length($read_sequence)-$number_of_S_letters_at_the_right_side_of_CIGAR,$number_of_S_letters_at_the_right_side_of_CIGAR));
				$average_length_of_extensions+=$number_of_S_letters_at_the_right_side_of_CIGAR;
			}
		}

	}
	close(SAMFILE);
	
	if($number_of_overhanging_reads_that_map_well_enough==0)
	{
			print LOGFILE "\nEnding, because there are no extensions\n";
			print HISTORY_OF_ALTERNATIVE_EXTENSIONS "\nElloreas has finished. See elloreas_logs.txt\n";
			open CONTIG_FROM_THE_FINAL_ITERATION, "> $output_folder/contig_from_the_final_iteration.fasta";
			print CONTIG_FROM_THE_FINAL_ITERATION ">$contig_title\n";
			print CONTIG_FROM_THE_FINAL_ITERATION "$contig_sequence\n";
			close(CONTIG_FROM_THE_FINAL_ITERATION);
			&make_the_final_coverage_plot_and_dot_plot();
			exit();
	}

	##########################################################################################################
	##########################################################################################################	

	#Now I, just for statistics, create the list of most frequent extensions
	
	#Теперь мы составляем хит-парад k-merов, вылезающих за край контига на k==$median_overhang_length
	%hash_number_of_such_extensions=(); #хэш с числом раз, которые встретился данный k-mer (ключи - собственно k-merы)
	foreach $this_read_extension_sequence(@array_extensions)
	{
		if(length($this_read_extension_sequence)>=$median_overhang_length)
		{
			$extending_kmer=substr($this_read_extension_sequence,0,$median_overhang_length);
			$hash_number_of_such_extensions{$extending_kmer}++;
		}
	}
	
	#теперь с помощью сортировки определяем, какой же k-mer здесь самый частый
	@array_for_sorted_hash_keys=sort {$hash_number_of_such_extensions{$b} <=> $hash_number_of_such_extensions{$a}} (keys %hash_number_of_such_extensions);
	
	#the number of alternative extensions
	$number_of_alternative_extensions=$#array_for_sorted_hash_keys+1;
=head
	print LOGFILE "There are $number_of_alternative_extensions alternative $median_overhang_length-mer extensions for the contig edge. Top is:\n";
	#теперь пишем хит-парад удлинений. Если удлинений больше 5, то пишем только самый частые 5.
	$number_of_extension=1;
	while(($number_of_extension<=5)&&($number_of_extension<=$number_of_alternative_extensions))
	{
		$k_mer_sequence=$array_for_sorted_hash_keys[$number_of_extension-1];
		print LOGFILE "top$number_of_extension is $k_mer_sequence which appears ".$hash_number_of_such_extensions{$k_mer_sequence}." times\n";
		$number_of_extension++;
	}
=cut	
	##########################################################################################################
	##########################################################################################################
	#Now I cluster k-mers (if the user provided $min_sequence_similarity<1) or take the most frequent k-mer (if the user provided $min_sequence_similarity=1)
	
	if($min_sequence_similarity<1)
	{
		#Now I create a FASTA-file with all overhangs which are at least $median_overhang_length long. I truncate them to the $median_overhang_length (thus producing k-mers). This file will be used for clustering with UCLUST.
		open FASTA_FILE_WITH_KMERS, "> $output_folder/Iteration$iteration_number/kmers_to_find_extensions.it$iteration_number.fasta";
		$this_overhang_number=0;
		foreach $this_read_extension_sequence(@array_extensions)
		{
			if(length($this_read_extension_sequence)>=$median_overhang_length)
			{
				$this_overhang_number++;
				$kmer=substr($this_read_extension_sequence,0,$median_overhang_length);
				print FASTA_FILE_WITH_KMERS ">kmer_$this_overhang_number\n$kmer\n";
			}
		}
		
		close(FASTA_FILE_WITH_KMERS);
		
		#cluster reads
		$similarity_for_clustering=$min_sequence_similarity*$min_sequence_similarity; #if we expect X% of similarity between a read and a reference then the expectation of similarity between two reads is approximately X*X. I'm not 100% sure, but I suppose it to be so
		system("usearch -cluster_fast $output_folder/Iteration$iteration_number/kmers_to_find_extensions.it$iteration_number.fasta -id $similarity_for_clustering -strand plus -threads $cores_to_use -sizeout -clusters $output_folder/Iteration$iteration_number/cluster");
		
		@array_of_consensus_sequences_without_Ns=(); #the array of consensus sequences where I have removed Ns (ambiguous bases. They usually arise in columns where some read has a sequencing error which resulted in an insertion. Also, they can arise in columns with so many deletional or point errors that EMBOSS cons cannot calculate the consensus). $array_of_consensus_sequences_without_Ns[0..4]="SEQUENCE" . Numbers of clusters from which these consensuses come are counted from 0 here, not from 1.
		
		#performing alignment and alignment consensus calculation for top-5 clusters (or less clusters, if there are less than 5). USEARCH outputs clusters in files with titles like cluster0, cluster1...
		
		print HISTORY_OF_ALTERNATIVE_EXTENSIONS "Alternative extensions most supported by reads:\n";
		print LOGFILE "Alternative extensions most supported by reads:\n";
		
		#now I determine the top-5 clusters represented by most reads.
		%hash_cluster_number_created_by_USEARCH_to_the_number_of_reads_forming_it=(); #$hash_cluster_number_created_by_USEARCH_to_the_number_of_reads_forming_it{17}=42 . ($hash_cluster_number_created_by_USEARCH_to_the_number_of_reads_forming_it{cluster_number_created_by_USEARCH}=how_many_reads_form_this_cluster) In this hash cluster numbers are counted from 0, not 1.
		opendir FOLDER_OF_THIS_ITERATION, "$output_folder/Iteration$iteration_number";
		@array_list_of_files_in_the_iteration_folder=readdir(FOLDER_OF_THIS_ITERATION);
		foreach $filename(@array_list_of_files_in_the_iteration_folder)
		{
			#if it is a FASTA file with kmers forming a cluster
			if($filename=~/cluster(\d+)$/)
			{
				$cluster_number_created_by_USEARCH=$1;
				open FILE_WITH_KMERS_FORMING_THE_CLUSTER, "< $output_folder/Iteration$iteration_number/cluster$cluster_number_created_by_USEARCH";
				$how_many_reads_support_this_extension=0;
				while(<FILE_WITH_KMERS_FORMING_THE_CLUSTER>)
				{
					if($_=~/^>/)
					{
						$how_many_reads_support_this_extension++;
					}
				}
				close(FILE_WITH_KMERS_FORMING_THE_CLUSTER);
				$hash_cluster_number_created_by_USEARCH_to_the_number_of_reads_forming_it{$cluster_number_created_by_USEARCH}=$how_many_reads_support_this_extension;
			}
		}
		closedir(FOLDER_OF_THIS_ITERATION);
		
		#creating an array with cluster numbers sorted in the order of decreasing number of reads forming them
		@array_of_cluster_numbers_created_by_USEARCH_sorted_by_decreasing_numbers_of_reads_forming_them=sort {$hash_cluster_number_created_by_USEARCH_to_the_number_of_reads_forming_it{$b} <=> $hash_cluster_number_created_by_USEARCH_to_the_number_of_reads_forming_it{$a}} (keys %hash_cluster_number_created_by_USEARCH_to_the_number_of_reads_forming_it);

		$USEARCH_number_of_the_largest_cluster=$array_of_cluster_numbers_created_by_USEARCH_sorted_by_decreasing_numbers_of_reads_forming_them[0];
		$size_of_the_largest_cluster=$hash_cluster_number_created_by_USEARCH_to_the_number_of_reads_forming_it{$USEARCH_number_of_the_largest_cluster};
		
		$maximum_number_of_cluster_in_the_top=0; #how many clusters are in the top. If USEARCH created more than 5 clusters, the top will contain only 5. If USEARCH created less than 5, the top will contain this number of clusters. $maximum_number_of_cluster_in_the_top is counted from 1, not 0.
		if(($#array_of_cluster_numbers_created_by_USEARCH_sorted_by_decreasing_numbers_of_reads_forming_them+1)>=5)
		{
			$maximum_number_of_cluster_in_the_top=5;
		}
		else
		{
			$maximum_number_of_cluster_in_the_top=$#array_of_cluster_numbers_created_by_USEARCH_sorted_by_decreasing_numbers_of_reads_forming_them+1;
		}
		
		foreach $number_of_cluster_in_the_top(1..$maximum_number_of_cluster_in_the_top) #1 is the cluster represented by most reads, 2 by less, 3 by even less, etc. IMPORTANT: this is not the same as a cluster number ($cluster_number_created_by_USEARCH). $cluster_number_created_by_USEARCH is the number assigned to the cluster by USEARCH, not the number in the top. $number_of_cluster_in_the_top is counted from 1, not 0.
		{
			
			$cluster_number_created_by_USEARCH=$array_of_cluster_numbers_created_by_USEARCH_sorted_by_decreasing_numbers_of_reads_forming_them[$number_of_cluster_in_the_top-1]; #counted from 0, not from 1
			
			#if this cluster is supported by only one read, then I don't need to perform a multiple alignment and to calculate a consensus. I just take the sequence of the respective k-mer as an alternative extension.
			if($hash_cluster_number_created_by_USEARCH_to_the_number_of_reads_forming_it{$cluster_number_created_by_USEARCH}==1)
			{
				open CLUSTER_FILE, "< $output_folder/Iteration$iteration_number/cluster$cluster_number_created_by_USEARCH";
				while(<CLUSTER_FILE>)
				{
					if($_!~/^>/)
					{
						$string=$_;
						chomp($string);
						$array_of_consensus_sequences_without_Ns[$cluster_number_created_by_USEARCH].=$string;
					}
				}
				close(CLUSTER_FILE);
			}
			else
			{
				$cluster_number_created_by_USEARCH=$array_of_cluster_numbers_created_by_USEARCH_sorted_by_decreasing_numbers_of_reads_forming_them[$number_of_cluster_in_the_top-1]; #counted from 0, not from 1
				
				system("mafft --localpair --maxiterate 1000 --thread $cores_to_use --nuc --quiet $output_folder/Iteration$iteration_number/cluster$cluster_number_created_by_USEARCH > $output_folder/Iteration$iteration_number/cluster$cluster_number_created_by_USEARCH"."_alignment.fasta");
				system("cons -sequence $output_folder/Iteration$iteration_number/cluster$cluster_number_created_by_USEARCH"."_alignment.fasta -datafile EDNAFULL -outseq $output_folder/Iteration$iteration_number/cluster$cluster_number_created_by_USEARCH"."_consensus.fasta");
				#loading sequence and removing Ns.
				open CONSENSUS_FILE, "< $output_folder/Iteration$iteration_number/cluster$cluster_number_created_by_USEARCH"."_consensus.fasta";
				while(<CONSENSUS_FILE>)
				{
					if($_!~/^>/)
					{
						$string=$_;
						$string=~s/[Nn]//g;
						$string=uc($string);
						chomp($string);
						$array_of_consensus_sequences_without_Ns[$cluster_number_created_by_USEARCH].=$string;
					}
				}
				close(CONSENSUS_FILE);
			}
			
			#removing the right 10% bases of the consensus, because the rightmost bases of the consensus often contain mistakes (see the comment "1)" to the version 1.15 in the file "my_comments_to_versions.txt").
			$length_of_the_consensus_before_chopping_the_right_10_percents=length($array_of_consensus_sequences_without_Ns[$cluster_number_created_by_USEARCH]);
			$array_of_consensus_sequences_without_Ns[$cluster_number_created_by_USEARCH]=substr($array_of_consensus_sequences_without_Ns[$cluster_number_created_by_USEARCH],0,int(0.9*$length_of_the_consensus_before_chopping_the_right_10_percents)+1); #"+1" because in an unusual case when the consensus was only 1 bp long, if I don't use "+1" I will elongate the contig by 0 bp, thus Elloreas will forever loop, continuing iteration with the same contig sequence, adding 0 bp each time.
			
			
			print HISTORY_OF_ALTERNATIVE_EXTENSIONS "Extension ".($number_of_cluster_in_the_top).") Supported by ".$hash_cluster_number_created_by_USEARCH_to_the_number_of_reads_forming_it{$cluster_number_created_by_USEARCH}." reads: $array_of_consensus_sequences_without_Ns[$cluster_number_created_by_USEARCH]\n";
			print LOGFILE "Extension ".($number_of_cluster_in_the_top).") Supported by ".$hash_cluster_number_created_by_USEARCH_to_the_number_of_reads_forming_it{$cluster_number_created_by_USEARCH}." reads\n";
		}


		$sequence_of_the_largest_cluster=$array_of_consensus_sequences_without_Ns[$USEARCH_number_of_the_largest_cluster];
		
		##########################################################################################################
		##########################################################################################################
		
		if($number_of_overhanging_reads_that_map_well_enough<$min_allowed_total_number_of_extending_kmers)
		{
			print LOGFILE "\n\nEnding, because there are less than $min_allowed_total_number_of_extending_kmers extensions\n";
			print HISTORY_OF_ALTERNATIVE_EXTENSIONS "\n\nElloreas has finished. See elloreas_logs.txt\n";
			open CONTIG_FROM_THE_FINAL_ITERATION, "> $output_folder/contig_from_the_final_iteration.fasta";
			print CONTIG_FROM_THE_FINAL_ITERATION ">$contig_title\n";
			print CONTIG_FROM_THE_FINAL_ITERATION "$contig_sequence\n";
			close(CONTIG_FROM_THE_FINAL_ITERATION);
			&make_the_final_coverage_plot_and_dot_plot();
			exit();
		}
		
		#if the number of reads forming the top cluster is less than $min_allowed_ratio_of_extensions*$number_of_overhanging_reads_that_map_well_enough, I'm stopping the script.
		if($size_of_the_largest_cluster<$min_allowed_ratio_of_extensions*$number_of_overhanging_reads_that_map_well_enough)
		{
			print LOGFILE "Ending, because of the insufficient ratio (".$size_of_the_largest_cluster."/".$number_of_overhanging_reads_that_map_well_enough.") of most frequent extension\n";
			print HISTORY_OF_ALTERNATIVE_EXTENSIONS "Elloreas has finished. See elloreas_logs.txt\n";
			open CONTIG_FROM_THE_FINAL_ITERATION, "> $output_folder/contig_from_the_final_iteration.fasta";
			print CONTIG_FROM_THE_FINAL_ITERATION ">$contig_title\n";
			print CONTIG_FROM_THE_FINAL_ITERATION "$contig_sequence\n";
			close(CONTIG_FROM_THE_FINAL_ITERATION);
			&make_the_final_coverage_plot_and_dot_plot();
			exit();
		}
		#if the number of reads forming the top cluster is less than $min_allowed_number_of_reads_supporting_the_most_popular_extension
		if($size_of_the_largest_cluster<$min_allowed_number_of_reads_supporting_the_most_popular_extension)
		{
			print LOGFILE "\n\nEnding, because of the insufficient number of reads ($size_of_the_largest_cluster) forming the largest cluster\n";
			print HISTORY_OF_ALTERNATIVE_EXTENSIONS "\n\nElloreas has finished. See elloreas_logs.txt\n";
			open CONTIG_FROM_THE_FINAL_ITERATION, "> $output_folder/contig_from_the_final_iteration.fasta";
			print CONTIG_FROM_THE_FINAL_ITERATION ">$contig_title\n";
			print CONTIG_FROM_THE_FINAL_ITERATION "$contig_sequence\n";
			close(CONTIG_FROM_THE_FINAL_ITERATION);
			&make_the_final_coverage_plot_and_dot_plot();
			exit();
		}
		
		#удлиняем последовательность контига
		$contig_sequence.=$sequence_of_the_largest_cluster;
	}
	else
	{
		#printing top frequent extensions
		print HISTORY_OF_ALTERNATIVE_EXTENSIONS "Alternative extensions most supported by reads:\n";
		print LOGFILE "Alternative extensions most supported by reads:\n";
		$number_of_extension=1;
		while(($number_of_extension<=5)&&($number_of_extension<=$number_of_alternative_extensions))
		{
			$k_mer_sequence=$array_for_sorted_hash_keys[$number_of_extension-1];
			print HISTORY_OF_ALTERNATIVE_EXTENSIONS "Extension $number_of_extension) Supported by $hash_number_of_such_extensions{$k_mer_sequence} reads: $k_mer_sequence\n";
			print LOGFILE "Extension $number_of_extension) Supported by $hash_number_of_such_extensions{$k_mer_sequence} reads\n";
			$number_of_extension++;
		}
		
		#if the frequency of the top k-mer was less than $min_allowed_ratio_of_extensions*$number_of_overhanging_reads_that_map_well_enough , I'm stopping the script.
		$top_kmer_sequence=$array_for_sorted_hash_keys[0];
		if($hash_number_of_such_extensions{$top_kmer_sequence}<$min_allowed_ratio_of_extensions*$number_of_overhanging_reads_that_map_well_enough)
		{
			print LOGFILE "\n\nEnding, because of the insufficient ratio (".$hash_number_of_such_extensions{$top_kmer_sequence}."/".$number_of_overhanging_reads_that_map_well_enough.") of most frequent extensions\n";
			print HISTORY_OF_ALTERNATIVE_EXTENSIONS "\n\nElloreas has finished. See elloreas_logs.txt\n";
			open CONTIG_FROM_THE_FINAL_ITERATION, "> $output_folder/contig_from_the_final_iteration.fasta";
			print CONTIG_FROM_THE_FINAL_ITERATION ">$contig_title\n";
			print CONTIG_FROM_THE_FINAL_ITERATION "$contig_sequence\n";
			close(CONTIG_FROM_THE_FINAL_ITERATION);
			&make_the_final_coverage_plot_and_dot_plot();
			exit();
		}
		#if the number of reads forming the top cluster is less than $min_allowed_number_of_reads_supporting_the_most_popular_extension
		if($hash_number_of_such_extensions{$top_kmer_sequence}<$min_allowed_number_of_reads_supporting_the_most_popular_extension)
		{
			print LOGFILE "\n\nEnding, because of the insufficient number of reads (".$hash_number_of_such_extensions{$top_kmer_sequence}.") forming the largest cluster\n";
			print HISTORY_OF_ALTERNATIVE_EXTENSIONS "\n\nElloreas has finished. See elloreas_logs.txt\n";
			open CONTIG_FROM_THE_FINAL_ITERATION, "> $output_folder/contig_from_the_final_iteration.fasta";
			print CONTIG_FROM_THE_FINAL_ITERATION ">$contig_title\n";
			print CONTIG_FROM_THE_FINAL_ITERATION "$contig_sequence\n";
			close(CONTIG_FROM_THE_FINAL_ITERATION);
			&make_the_final_coverage_plot_and_dot_plot();
			exit();
		}
		
		if($number_of_overhanging_reads_that_map_well_enough<$min_allowed_total_number_of_extending_kmers)
		{
			print LOGFILE "\n\nEnding, because there are less than $min_allowed_total_number_of_extending_kmers extensions\n";
			print HISTORY_OF_ALTERNATIVE_EXTENSIONS "\n\nElloreas has finished. See elloreas_logs.txt\n";
			open CONTIG_FROM_THE_FINAL_ITERATION, "> $output_folder/contig_from_the_final_iteration.fasta";
			print CONTIG_FROM_THE_FINAL_ITERATION ">$contig_title\n";
			print CONTIG_FROM_THE_FINAL_ITERATION "$contig_sequence\n";
			close(CONTIG_FROM_THE_FINAL_ITERATION);
			&make_the_final_coverage_plot_and_dot_plot();
			exit();
		}
		
		#удлиняем последовательность контига
		$contig_sequence.=$array_for_sorted_hash_keys[0];
	}
	
	$new_contig_filename=$starter_file_path;
	$new_contig_filename=~s/^(.+)\..*?$/\1.it$iteration_number.fasta/;
	#теперь создаём файл с удлинённым контигом
	open CURRENT_CONTIG_OUTFILE, "> $output_folder/Iteration$iteration_number/$new_contig_filename";
	print CURRENT_CONTIG_OUTFILE ">$contig_title\n$contig_sequence\n";
	close(CURRENT_CONTIG_OUTFILE);
	
	#now, if the user asked for this, I remove read mappings (sam, bam and bam.bai files)
	if($should_elloreas_delete_sam_and_bam_files=~/true/)
	{
		system("rm -f $output_folder/Iteration$iteration_number/mapping.it$iteration_number.sam");
		system("rm -f $output_folder/Iteration$iteration_number/mapping.it$iteration_number.bam");
		system("rm -f $output_folder/Iteration$iteration_number/mapping.it$iteration_number.sorted.bam");
		system("rm -f $output_folder/Iteration$iteration_number/mapping.it$iteration_number.sorted.bam.bai");
	}
	
	$current_contig_filename=$new_contig_filename;
	
	#now check with BLAST whether Elloreas has a terminal sequence which causes eternal looping.
	system("blastn -task megablast -query $output_folder/Iteration$iteration_number/$current_contig_filename -subject $output_folder/Iteration$iteration_number/$current_contig_filename -out $output_folder/Iteration$iteration_number/megablast_results.it$iteration_number.txt -outfmt \"6 qseqid sseqid evalue pident qstart qend sstart send length sstrand\" -num_threads $cores_to_use -max_target_seqs 1 -max_hsps 1000000 -evalue 1e-10");
	open BLAST_RESULTS, "< $output_folder/Iteration$iteration_number/megablast_results.it$iteration_number.txt";
	@array_file_with_blast_results=<BLAST_RESULTS>;
=head
Paenibacillus_sp_RUD330__starting_from_dnaA	Paenibacillus_sp_RUD330__starting_from_dnaA	0.0	100.000	1	31825	1	31825	31825	plus
Paenibacillus_sp_RUD330__starting_from_dnaA	Paenibacillus_sp_RUD330__starting_from_dnaA	0.0	81.068	26531	31721	21236	26451	5467	plus
=cut
	#if there is only one BLAST result, then there are no repeats at all, because the first line ($array_file_with_blast_results[0]) is always the alignment of the entire contig to itself.
	if($array_file_with_blast_results[1])
	{
		$line_number=1;
		#print LOGFILE, "Blast results string for the repeat with the highest e-value: $array_file_with_blast_results[1] (even if the looping repeat exists, it is not always the looping repeat";
		while($array_file_with_blast_results[$line_number])
		{
			#$array_string_split[9]=~/plus/ is checked because inverted repeats do not cause loops.
			#also checking that the repeat is perfect (100% sequence similarity)
			@array_string_split=split(/\t/,$array_file_with_blast_results[$line_number]);
			if((($array_string_split[5]==length($contig_sequence))||($array_string_split[7]==length($contig_sequence)))&&($array_string_split[9]=~/plus/)&&($array_string_split[8]>=$min_aligned_length)&&($array_string_split[3]==100))
			{
				print LOGFILE "\n\nEnding, because Elloreas has entered an eternal loop. This means one of the following two things:\n";
				print LOGFILE "1) You have assembled a circular genome.\n";
				print LOGFILE "2) There is a large repeat with one repeat unit somewhere within the contig and the other repeat unit at the right (3') edge of the contig.\n";
				print LOGFILE "To differentiate between these two cases you can look at the dot plot created by Elloreas (dot_plot_for_the_contig_from_the_final_iteration.jpeg). If you see one repeat unit at the left (5') edge of the contig and the other at the right edge (3'), this means that the contig has circularized (\"1)\"). Basically, during the elongation process Elloreas has started creating on its right edge the same sequence as is on its left edge. If this is a circular genome, you can now remove one of these repeat instances.\n";
				print LOGFILE "If you see that the left repeat unit is not at the left edge of the contig, this means you have encountered a long repeat present in your genome (\"2)\"). Letting Elloreas work further will make him fall into an eternal loop, forever looping between these two repeat units. You may want to look in the file history_of_alternative_extensions.txt to find a highly supported alternative extension on one of recent iterations of Elloreas, and then choose this alternative extension. See the part of the Frequently Asked Questions called \"There was a fork at some iteration. Elloreas chose one of alternative extensions, but I want to try another one. What should I do?\"\n";
				print HISTORY_OF_ALTERNATIVE_EXTENSIONS "Elloreas has finished. See elloreas_logs.txt\n";
				open CONTIG_FROM_THE_FINAL_ITERATION, "> $output_folder/contig_from_the_final_iteration.fasta";
				print CONTIG_FROM_THE_FINAL_ITERATION ">$contig_title\n";
				print CONTIG_FROM_THE_FINAL_ITERATION "$contig_sequence\n";
				close(CONTIG_FROM_THE_FINAL_ITERATION);
				&make_the_final_coverage_plot_and_dot_plot();
				exit();
			}
			else
			{
				#print LOGFILE "Not stopping Elloreas, BLAST line is $array_file_with_blast_results[$line_number] , values are ".$array_string_split[5]." ".length($contig_sequence)." ".$array_string_split[7]." ".length($contig_sequence)." ".$array_string_split[9]." ".$array_string_split[8]." ".$min_aligned_length."\n";
			}
			$line_number++;
		}
		
	}
	else
	{
		#print LOGFILE, "Self-BLAST suggests that there are no repeats at all\n";
	}
	close(BLAST_RESULTS);
	
}


print LOGILE "\nEnding, because Elloreas has reached the final iteration (number $max_number_of_cycles) according to the value provided by the user with the \"--maximum_number_of_iterations\" option.\n";


#this subroutine takes the sequence of the final contig produced by Elloreas and makes a coverage plot and a dot plot. Temporary files are deleted.
sub make_the_final_coverage_plot_and_dot_plot()
{
		system("mkdir $output_folder/Temporary");
		system("minimap2 -x $sequencing_technology -a --secondary no --sam-hit-only -o $output_folder/Temporary/mapping.sam -t $cores_to_use $output_folder/contig_from_the_final_iteration.fasta $reads_file_path");
		system("samtools view -Sbh $output_folder/Temporary/mapping.sam >$output_folder/Temporary/mapping.bam");
		system("samtools sort $output_folder/Temporary/mapping.bam >$output_folder/Temporary/mapping.sorted.bam");
		system("samtools index $output_folder/Temporary/mapping.sorted.bam");
		system("bedtools genomecov -d -ibam $output_folder/Temporary/mapping.sorted.bam > $output_folder/Temporary/coverage_list.txt");
		system("Rscript $path_to_the_folder_where_Elloreas_is_located/plot_coverage_for_Elloreas.r $output_folder/Temporary/coverage_list.txt $output_folder/coverage_plot_for_the_contig_from_the_final_iteration.jpeg");
		#making the dotplot
		system("lastz $output_folder/contig_from_the_final_iteration.fasta $output_folder/contig_from_the_final_iteration.fasta --notransition --strand=both --step=20 --nogapped --format=rdotplot > $output_folder/Temporary/data_to_build_the_dotplot.txt");
		system("Rscript $path_to_the_folder_where_Elloreas_is_located/make_dotplot_for_Elloreas.r $output_folder/Temporary/data_to_build_the_dotplot.txt $output_folder/dotplot_for_the_contig_from_the_final_iteration.jpeg");
		#removing the temporary folder
		system("rm -rf $output_folder/Temporary");
		#if the user asked for this, remove data from all iterations which haven't been deleted before (basically, its only the very last iteration)
		if($should_elloreas_delete_folders_of_all_iterations=~/true/)
		{
			system("rm -rf ./$output_folder/Iteration*");
		}
}











################################################################################################################
################################################################################################################
#Subroutines
################################################################################################################

#the subroutine to take parameters provided by the user in a command line
sub take_parameters_from_the_command_line()
{
	#assigning default values to parameters
	$starter_file_path=""; #one of the two parameters that have no default values
	$reads_file_path=""; #one of the two parameters that have no default values
	$min_sequence_similarity="80%";
	$min_aligned_length=5000;
	$contig_edge_size_to_use_for_alignment=10000;
	$min_allowed_total_number_of_extending_kmers=4;
	$min_allowed_ratio_of_extensions="1%";
	$max_number_of_cycles=1000000;
	$cores_to_use=10;
	$number_of_bases_to_trim_from_the_starter_edges_in_the_very_beginning_of_the_assembly=100;
	$output_folder="Elloreas_output";
	$min_allowed_number_of_reads_supporting_the_most_popular_extension=4;
	$should_elloreas_draw_coverage_plots="false";
	$sequencing_technology="oxford_nanopore";
	$should_elloreas_delete_sam_and_bam_files="true";
	$should_elloreas_delete_folders_of_all_iterations="false";

	#$hash_list_of_allowed_keys_to_the_word_yes{"--reads"}="yes". Only the keys present in this hash are accepted from the user. If the user provides a key which is absent from this hash, Elloreas prints an error and stops. This hash has no values other than "yes".
	%hash_list_of_allowed_keys_to_the_word_yes=(
	"--reads"=>"yes",
	"--starter"=>"yes",
	"--sequencing_technology"=>"yes",
	"--minimum_length_of_mapped_read_part"=>"yes",
	"--number_of_CPU_threads"=>"yes",
	"--output_folder"=>"yes",
	"--minimum_read_similarity"=>"yes",
	"--contig_edge_size_to_use_for_mapping"=>"yes",
	"--minimum_total_number_of_extending_reads"=>"yes",
	"--minimum_percent_of_reads_in_the_most_popular_extension"=>"yes",
	"--minimum_number_of_reads_in_the_most_popular_extension"=>"yes",
	"--maximum_number_of_iterations"=>"yes",
	"--number_of_bases_to_trim_from_the_starter_edges"=>"yes",
	"--draw_coverage_plots"=>"yes",
	"--delete_sam_and_bam"=>"yes",
	"--delete_folders_of_iterations"=>"yes",
	"--help"=>"yes",
	"--version"=>"yes");
	
	#this is the content of Elloreas's help
$the_content_of_Elloreas_help=<<HERE_DOCUMENT;
=====================================================================
=====================================================================
Elloreas, a genome assembler. For the complete manual, read the manual.pdf file which is located in the same folder as the executable of Elloreas. You can download the latest version of Elloreas from https://github.com/shelkmike/Elloreas . 
=====================================================================
=====================================================================

Only two options of Elloreas are mandatory: --reads and --starter.

Main options:
1) --reads		Path to the file with sequencing reads. They can be in FASTQ or FASTA, gzipped or not. Example: /mnt/some_folder/reads.fastq

2) --starter		Path to the file with the contig that you want to elongate (alternatively, you can put there a sequence of a random read instead of a contig). This should be a FASTA file containing only one sequence. Example: /mnt/some_folder/starter.fasta

3) --sequencing_technology		Which sequencing technology these reads come from. Three variants are allowed: typical_pacbio, hifi_pacbio, oxford_nanopore. The default value is oxford_nanopore. The option \"hifi_pacbio\" if for reads produced with the PacBio CCS technology.

4) --minimum_length_of_mapped_read_part		How many bases from a read should map to the contig, so Elloreas considers this read. The default value is 5000. I recommend to set this value approximately to 75% of the median length of input reads.

5) --number_of_CPU_threads		How many CPU threads to use. The default value is 10.

=====================================================================
=====================================================================

Less important options (usually they shouldn't be changed):

6) --output_folder		The folder where the results should be printed. The default value is ./Elloreas_output

7) --minimum_read_similarity		Minimum percent of sequence similarity between the sequence of a read and the sequence of the contig. Only reads with at least this similarity will be considered by Elloreas for extension. Examples: 80% or 99.7%. Default value: 80% if the user specifies typical_pacbio or oxford_nanopore as the sequencing technology (see above) or 98% if the user specifies hifi_pacbio as the sequencing technology. If the user does not specify the --sequencing_technology option at all, this value is set to 80%.

8) --contig_edge_size_to_use_for_mapping		How many bases on each iteration to take from the contig's edge for mapping. I recommend to set this approximately to the median read length multiplied by 1.5. Example: 10000. The default value is 10000 or contig's length, whichever is smaller.

9) --minimum_total_number_of_extending_reads		Elloreas will stop if there are less than this number of reads that overhang the contig's edge. Example: 4. The default value is 4.

10) --minimum_percent_of_reads_in_the_most_popular_extension		Elloreas will stop if the most popular extension is supported by less than this percent of reads. Examples: 40% or 1%. The default value is 1%. Because of this low default value, Elloreas by default will not stop even at forks with a huge amount of alternative extensions.

11) --minimum_number_of_reads_in_the_most_popular_extension		Similar to the previous option, but the number of reads instead of the percent of reads. Elloreas will stop if the most popular extension is supported by less than this number of reads. Examples: 10 or 4. The default value is 4.

12) --maximum_number_of_iterations		How many iterations of elongation Elloreas performs before it stops. The default value is 1000000. Basically, such a high value removes the criterion on the maximum number of iterations.

13) --number_of_bases_to_trim_from_the_starter_edges		How many bases to trim from both edges of the starter at the very beginning of the work of Elloreas. Trimming some bases is useful because there may be assembly errors at the edge of the starter. The default value is 100.

14) --draw_coverage_plots		Should coverage plots be drawn for contig edges at each iteration? Possible values: true or false. The default value is false. Making coverage plots slightly increases the computation time of Elloreas.

15) --delete_sam_and_bam		Should sam and bam files produced during each iteration be kept? Possible values: true or false. The default value is true. Keeping these files increases the disk space required for Elloreas.

16) --delete_folders_of_iterations		Should Elloreas delete folders where it stores data corresponding to different iterations? Possible values: true or false. The default value is false. Switching this to true will safe a bit of disk space.

17) --help		print this help

18) --version		print the version of Elloreas

=====================================================================
=====================================================================

A simple example:
perl elloreas.pl --reads reads.fastq --starter starter.fasta

A not so simple example:
perl elloreas.pl --reads reads.fastq --starter starter.fasta --sequencing_technology oxford_nanopore --minimum_length_of_aligned_read_part 8000

A more elaborate example:
perl /home/user/Work/Tools/Elloreas/elloreas.pl --reads /home/user/data/Reads/RLT1742.fastq --starter /home/user/data/Assembly/contig_1_from_Spades_assembly.fasta --output_folder /home/user/data/Assembly/Elloreas/Output --sequencing_technology typical_pacbio --minimum_length_of_aligned_read_part 8000 --number_of_CPU_threads 20 --delete_sam_and_bam false --minimum_read_similarity 90%
		
=====================================================================
=====================================================================
HERE_DOCUMENT

	$list_of_errors_in_the_input_command__to_print="Unfortunately, there are errors in the input:\n"; #if a user made errors in the input command, they all are accumulated in this variable and then printed (after all such errors are found).
	$number_of_the_current_error=0; #all errors are enumerated starting from 1
	#number of elements in the command that the user provided
	$number_of_elements_in_the_command=$#ARGV+1;
	
	#if the user has provided no options, I just print help and stop Elloreas;
	if($number_of_elements_in_the_command==0)
	{
		print $the_content_of_Elloreas_help;
		exit();
	}
	
	#checking if --help or --version are present among the elements. They are the only two keys that don't have values.
	foreach $element (@ARGV)
	{
		#if that element was the last element in the command line, it may contain the newline character. I should trim it
		chomp($element);
		
		if($element=~/--help/)
		{
			print $the_content_of_Elloreas_help;
			exit();
		}
		if($element=~/--version/)
		{
			print "The version of Elloreas is $Elloreas_version\n";
			exit();
		}
	}

	$current_element_number=1;
	%hash_input_key_to_input_value=(); #like {"--starter"}="/mnt/some_folder/starter.fasta"
	while($current_element_number<=$number_of_elements_in_the_command)
	{
		$key=$ARGV[$current_element_number-1];
		$value=$ARGV[$current_element_number];

		$hash_input_key_to_input_value{$key}=$value;	
		$current_element_number+=2;
		
		#if the user has provided a value (not a key) starting with "--", this means he made a serious mistake (presumably involving a shift in key names) which probably makes impossible parsing the rest of commands and finding errors among them. This is why Elloreas should stop immediately after finding such a mistake.
		if($value=~/^\-\-/)
		{
			print "Error: the value for the option $key is \"$value\". You've probably made a mistake in the command line.\n";
			exit();
		}
		
		#checking that Elloreas knows this key
		if($hash_list_of_allowed_keys_to_the_word_yes{$key}=~/^$/)
		{
			$number_of_the_current_error++;
			$list_of_errors_in_the_input_command__to_print.=$number_of_the_current_error.") You've entered a key which Elloreas doesn't know: \"$key\"\n";
		}
	}


	if($hash_input_key_to_input_value{"--starter"}!~/^$/)
	{
		$starter_file_path=$hash_input_key_to_input_value{"--starter"};
		
		#if this file does not exist
		if(! -e $starter_file_path)
		{
			$number_of_the_current_error++;
			$list_of_errors_in_the_input_command__to_print.=$number_of_the_current_error.") The file with the starter that you have specified ($starter_file_path) does not exist \n";
		}
	}	
	if($hash_input_key_to_input_value{"--reads"}!~/^$/)
	{
		$reads_file_path=$hash_input_key_to_input_value{"--reads"};

		#if this file does not exist
		if(! -e $reads_file_path)
		{
			$number_of_the_current_error++;
			$list_of_errors_in_the_input_command__to_print.=$number_of_the_current_error.") The file with reads that you have specified ($reads_file_path) does not exist\n";
		}
	}	
	if($hash_input_key_to_input_value{"--minimum_read_similarity"}!~/^$/)
	{
		$min_sequence_similarity=$hash_input_key_to_input_value{"--minimum_read_similarity"};
		if(($min_sequence_similarity!~/\d+\.\d+\%/)&&($min_sequence_similarity!~/\d+\%/))
		{
			$number_of_the_current_error++;
			$list_of_errors_in_the_input_command__to_print.=$number_of_the_current_error.") You have provided an improperly formatted value with the \"--minimum_read_similarity\" key. It should be formatted as 80%, while your is \"$min_sequence_similarity\"\n";		
		}
	}
	if($hash_input_key_to_input_value{"--minimum_length_of_mapped_read_part"}!~/^$/)
	{
		$min_aligned_length=$hash_input_key_to_input_value{"--minimum_length_of_mapped_read_part"};
		if($min_aligned_length!~/\d+/)
		{
			$number_of_the_current_error++;
			$list_of_errors_in_the_input_command__to_print.=$number_of_the_current_error.") You have provided an improperly formatted value with the \"--minimum_length_of_mapped_read_part\" key. It should be formatted as 2000, while your is \"$min_aligned_length\"\n";		
		}
	}
	if($hash_input_key_to_input_value{"--contig_edge_size_to_use_for_mapping"}!~/^$/)
	{
		$contig_edge_size_to_use_for_alignment=$hash_input_key_to_input_value{"--contig_edge_size_to_use_for_mapping"};
		if($contig_edge_size_to_use_for_alignment!~/\d+/)
		{
			$number_of_the_current_error++;
			$list_of_errors_in_the_input_command__to_print.=$number_of_the_current_error.") You have provided an improperly formatted value with the \"--contig_edge_size_to_use_for_mapping\" key. It should be formatted as 20000, while your is \"$contig_edge_size_to_use_for_alignment\"\n";
		}
	}
	if($hash_input_key_to_input_value{"--minimum_total_number_of_extending_reads"}!~/^$/)
	{
		$min_allowed_total_number_of_extending_kmers=$hash_input_key_to_input_value{"--minimum_total_number_of_extending_reads"};
		if($min_allowed_total_number_of_extending_kmers!~/\d+/)
		{
			$number_of_the_current_error++;
			$list_of_errors_in_the_input_command__to_print.=$number_of_the_current_error.") You have provided an improperly formatted value with the \"--minimum_total_number_of_extending_reads\" key. It should be formatted as 5, while your is \"$min_allowed_total_number_of_extending_kmers\"\n";
		}
	}
	if($hash_input_key_to_input_value{"--minimum_percent_of_reads_in_the_most_popular_extension"}!~/^$/)
	{
		$min_allowed_ratio_of_extensions=$hash_input_key_to_input_value{"--minimum_percent_of_reads_in_the_most_popular_extension"};
		if(($min_allowed_ratio_of_extensions!~/\d+\.\d+\%/)&&($min_allowed_ratio_of_extensions!~/\d+\%/))
		{
			$number_of_the_current_error++;
			$list_of_errors_in_the_input_command__to_print.=$number_of_the_current_error.") You have provided an improperly formatted value with the \"--minimum_percent_of_reads_in_the_most_popular_extension\" key. It should be formatted as 10%, while your is \"$min_allowed_ratio_of_extensions\"\n";
		}
	}
	if($hash_input_key_to_input_value{"--maximum_number_of_iterations"}!~/^$/)
	{
		$max_number_of_cycles=$hash_input_key_to_input_value{"--maximum_number_of_iterations"};
		if($max_number_of_cycles!~/\d+/)
		{
			$number_of_the_current_error++;
			$list_of_errors_in_the_input_command__to_print.=$number_of_the_current_error.") You have provided an improperly formatted value with the \"--maximum_number_of_iterations\" key. It should be formatted as 1000, while your is \"$max_number_of_cycles\"\n";
		}
	}
	if($hash_input_key_to_input_value{"--number_of_CPU_threads"}!~/^$/)
	{
		$cores_to_use=$hash_input_key_to_input_value{"--number_of_CPU_threads"};
		if($cores_to_use!~/\d+/)
		{
			$number_of_the_current_error++;
			$list_of_errors_in_the_input_command__to_print.=$number_of_the_current_error.") You have provided an improperly formatted value with the \"--number_of_CPU_threads\" key. It should be formatted as 1000, while your is \"$cores_to_use\"\n";
		}
	}
	if($hash_input_key_to_input_value{"--number_of_bases_to_trim_from_the_starter_edges"}!~/^$/)
	{
		$number_of_bases_to_trim_from_the_starter_edges_in_the_very_beginning_of_the_assembly=$hash_input_key_to_input_value{"--number_of_bases_to_trim_from_the_starter_edges"};
		if($number_of_bases_to_trim_from_the_starter_edges_in_the_very_beginning_of_the_assembly!~/\d+/)
		{
			$number_of_the_current_error++;
			$list_of_errors_in_the_input_command__to_print.=$number_of_the_current_error.") You have provided an improperly formatted value with the \"--number_of_bases_to_trim_from_the_starter_edges\" key. It should be formatted as 5000, while your is \"$number_of_bases_to_trim_from_the_starter_edges_in_the_very_beginning_of_the_assembly\"\n";
		}
	}
	if($hash_input_key_to_input_value{"--output_folder"}!~/^$/)
	{
		$output_folder=$hash_input_key_to_input_value{"--output_folder"};
	}
	if($hash_input_key_to_input_value{"--minimum_number_of_reads_in_the_most_popular_extension"}!~/^$/)
	{
		$min_allowed_number_of_reads_supporting_the_most_popular_extension=$hash_input_key_to_input_value{"--minimum_number_of_reads_in_the_most_popular_extension"};
		if($min_allowed_number_of_reads_supporting_the_most_popular_extension!~/\d+/)
		{
			$number_of_the_current_error++;
			$list_of_errors_in_the_input_command__to_print.=$number_of_the_current_error.") You have provided an improperly formatted value with the \"--minimum_number_of_reads_in_the_most_popular_extension\" key. It should be formatted as 2, while your is \"$min_allowed_number_of_reads_supporting_the_most_popular_extension\"\n";
		}
	}
	if($hash_input_key_to_input_value{"--draw_coverage_plots"}!~/^$/)
	{
		$should_elloreas_draw_coverage_plots=$hash_input_key_to_input_value{"--draw_coverage_plots"};
		if($should_elloreas_draw_coverage_plots!~/(true|false)/)
		{
			$number_of_the_current_error++;
			$list_of_errors_in_the_input_command__to_print.=$number_of_the_current_error.") You have provided an improper value with the \"--draw_coverage_plots\" key. It should be either true or false, while your is \"$should_elloreas_draw_coverage_plots\"\n";
		}
	}
	if($hash_input_key_to_input_value{"--sequencing_technology"}!~/^$/)
	{
		$sequencing_technology=$hash_input_key_to_input_value{"--sequencing_technology"};
		if($sequencing_technology!~/(typical_pacbio|hifi_pacbio|oxford_nanopore)/)
		{
			$number_of_the_current_error++;
			$list_of_errors_in_the_input_command__to_print.=$number_of_the_current_error.") You have provided an improper value with the \"--sequencing_technology\" key. It should be either typical_pacbio or hifi_pacbio or oxford_nanopore, while your is \"$sequencing_technology\"\n";
		}
	}
	if($hash_input_key_to_input_value{"--delete_sam_and_bam"}!~/^$/)
	{
		$should_elloreas_delete_sam_and_bam_files=$hash_input_key_to_input_value{"--delete_sam_and_bam"};
		if($should_elloreas_delete_sam_and_bam_files!~/(true|false)/)
		{
			$number_of_the_current_error++;
			$list_of_errors_in_the_input_command__to_print.=$number_of_the_current_error.") You have provided an improper value with the \"--delete_sam_and_bam\" key. It should be either true or false, while your is \"$should_elloreas_delete_sam_and_bam_files\"\n";
		}
	}
	if($hash_input_key_to_input_value{"--delete_folders_of_iterations"}!~/^$/)
	{
		$should_elloreas_delete_folders_of_all_iterations=$hash_input_key_to_input_value{"--delete_folders_of_iterations"};
		if($should_elloreas_delete_folders_of_all_iterations!~/(true|false)/)
		{
			$number_of_the_current_error++;
			$list_of_errors_in_the_input_command__to_print.=$number_of_the_current_error.") You have provided an improper value with the \"--delete_folders_of_iterations\" key. It should be either true or false, while your is \"$should_elloreas_delete_folders_of_all_iterations\"\n";
		}
	}
		
		

	#checking the mandatory parameters (the path to reads and the path to a starter)
	if($reads_file_path=~/^$/)
	{
		$number_of_the_current_error++;
		$list_of_errors_in_the_input_command__to_print.=$number_of_the_current_error.") Path to the file with reads is not specified. To do this, use a command like --reads /mnt/some_folder/reads.fastq \n";
	}
	if($starter_file_path=~/^$/)
	{
		$number_of_the_current_error++;
		$list_of_errors_in_the_input_command__to_print.=$number_of_the_current_error.") Path to the file with a starting sequence (\"starter\") is not specified. To do this, use a command like --starter /mnt/some_folder/starter.fasta \n";
	}

	#converting values that user provided as percents into decimals (like 80.2% to 0.802)
	$min_sequence_similarity=~s/\%$//;
	$min_sequence_similarity=$min_sequence_similarity/100;
	
	$min_allowed_ratio_of_extensions=~s/\%$//;
	$min_allowed_ratio_of_extensions=$min_allowed_ratio_of_extensions/100;
	

	#if the user did specify --sequencing_technology but didn't specify --minimum_read_similarity , I set --minimum_read_similarity to 80%.
	if(($hash_input_key_to_input_value{"--sequencing_technology"}!~/^$/)&&($hash_input_key_to_input_value{"--minimum_read_similarity"}=~/^$/))
	{
		if(($sequencing_technology=~/^typical_pacbio$/)||($sequencing_technology=~/^oxford_nanopore$/))
		{
			$min_sequence_similarity=0.8;
		}
		if($sequencing_technology=~/^hifi_pacbio$/)
		{
			$min_sequence_similarity=0.98;
		}
	}

	#converting the name of the sequencing technology to the form that can be understood by minimap2
	if($sequencing_technology=~/^typical_pacbio$/)
	{
		$sequencing_technology="map-pb";
	}
	if($sequencing_technology=~/^oxford_nanopore$/)
	{
		$sequencing_technology="map-ont";
	}
	if($sequencing_technology=~/^hifi_pacbio$/)
	{
		$sequencing_technology="asm20";
	}

	#If there exist errors in the input, Elloreas prints them and exits. It doesn't print them to "elloreas_logs.txt", because the path to the output folder is one of the parameters. Instead, Elloreas prints errors to the standard output.
	if($number_of_the_current_error>=1)
	{
		print "$list_of_errors_in_the_input_command__to_print\n";
		exit();
	}
}

#this subroutine checks whether programs required by Elloreas exist in $PATH and whether there are scripts plot_coverage_for_Elloreas.r and make_dotplot_for_Elloreas.r in the same folder where elloreas.pl lies.
sub check_whether_required_parameters_and_scripts_exist()
{
	$are_there_missing_programs="no"; #"no" if not, "yes" if there are missing programs
	print LOGFILE "Checking whether all programs required by Elloreas are present.\n";
	if(`which perl`=~/perl$/)
	{
		print LOGFILE "perl is present\n";
	}
	else
	{
		print LOGFILE "perl is absent\n";
		$are_there_missing_programs="yes";
	}

	if(`which Rscript`=~/Rscript$/)
	{
		print LOGFILE "Rscript is present\n";
	}
	else
	{
		print LOGFILE "Rscript is absent\n";
		$are_there_missing_programs="yes";
	}

	if(`which minimap2`=~/minimap2$/)
	{
		print LOGFILE "minimap2 is present\n";
	}
	else
	{
		print LOGFILE "minimap2 is absent\n";
		$are_there_missing_programs="yes";
	}

	if(`which samtools`=~/samtools$/)
	{
		print LOGFILE "samtools is present\n";
	}
	else
	{
		print LOGFILE "samtools is absent\n";
		$are_there_missing_programs="yes";
	}

	if(`which blastn`=~/blastn$/)
	{
		print LOGFILE "blastn is present\n";
	}
	else
	{
		print LOGFILE "blastn is absent\n";
		$are_there_missing_programs="yes";
	}

	if(`which bedtools`=~/bedtools$/)
	{
		print LOGFILE "bedtools is present\n";
	}
	else
	{
		print LOGFILE "bedtools is absent\n";
		$are_there_missing_programs="yes";
	}

	if(`which usearch`=~/usearch$/)
	{
		print LOGFILE "usearch is present\n";
	}
	else
	{
		print LOGFILE "usearch was not found. Sometimes executable files of USEARCH have complicated names, like \"usearch11.0.667_i86linux32\". Please, rename the executable file of USEARCH to just \"usearch\"\n";
		$are_there_missing_programs="yes";
	}

	if(`which mafft`=~/mafft$/)
	{
		print LOGFILE "mafft is present\n";
	}
	else
	{
		print LOGFILE "mafft was not found. Please, check that the path to the file \"mafft\" (not \"mafft.bat\" or something like this) is in the environment variable \$PATH\n";
		$are_there_missing_programs="yes";
	}

	if(`which cons`=~/cons$/)
	{
		print LOGFILE "EMBOSS is present\n";
	}
	else
	{
		print LOGFILE "EMBOSS was not found. Please, check that the path to its folder with binaries is present in the environment variable \$PATH. Namely, Elloreas requires the program \"cons\" that lies in this folder.\n";
		$are_there_missing_programs="yes";
	}
	
	if(`which lastz`=~/lastz$/)
	{
		print LOGFILE "lastz is present\n";
	}
	else
	{
		print LOGFILE "lastz is absent\n";
		$are_there_missing_programs="yes";
	}	
	
	#now checking the presence of scripts plot_coverage_for_Elloreas.r and make_dotplot_for_Elloreas.r which should lie in the same folder as elloreas.pl
	$presumed_path_to__plot_coverage_for_Elloreas=$path_to_the_folder_where_Elloreas_is_located."/plot_coverage_for_Elloreas.r";
	$presumed_path_to__plot_coverage_for_Elloreas=~s/\/\//\//g; #replacing // in the path by / . This doesn't change the accessibility of this script, but looks better.
	if(-e $presumed_path_to__plot_coverage_for_Elloreas)
	{
		
		print LOGFILE "$presumed_path_to__plot_coverage_for_Elloreas is present\n";
	}
	else
	{
		print LOGFILE "plot_coverage_for_Elloreas.r is absent. It is provided with Elloreas and should lie in the same folder where elloreas.pl is located ($path_to_the_folder_where_Elloreas_is_located). If you cannot find it, you can download it with Elloreas from https://github.com/shelkmike/Elloreas/releases.\n";
		$are_there_missing_programs="yes";
	}
	
	$presumed_path_to__make_dotplot_for_Elloreas=$path_to_the_folder_where_Elloreas_is_located."/make_dotplot_for_Elloreas.r";
	$presumed_path_to__make_dotplot_for_Elloreas=~s/\/\//\//g; #replacing // in the path by / . This doesn't change the accessibility of this script, but looks better.
	if(-e $presumed_path_to__make_dotplot_for_Elloreas)
	{
		print LOGFILE "$presumed_path_to__make_dotplot_for_Elloreas is present\n";
	}
	else
	{
		print LOGFILE "make_dotplot_for_Elloreas.r is absent. It is provided with Elloreas and should lie in the same folder where elloreas.pl is located ($path_to_the_folder_where_Elloreas_is_located). If you cannot find it, you can download it with Elloreas from https://github.com/shelkmike/Elloreas/releases.\n";
		$are_there_missing_programs="yes";
	}	
	
	if($are_there_missing_programs=~/^no$/)
	{
		print LOGFILE "Nice, all programs are here!\n\n";
	}
	else
	{
		print LOGFILE "\nUnfortunately, there are missing programs. They should be installed and added to the \$PATH environment variable.\n";
		exit();
	}
}

#this subroutine prints a list of parameters provided by the user via the command line. It uses the hash $hash_input_key_to_input_value{$key} created by the subroutine &take_parameters_from_the_command_line
sub print_parameters_provided_by_the_user_via_the_command_line()
{
	print LOGFILE "You have provided the following options:\n";
	$number_of_entered_option=0; 
	foreach $key(keys(%hash_input_key_to_input_value))
	{
		$number_of_entered_option++;
		print LOGFILE "$number_of_entered_option) $key ".$hash_input_key_to_input_value{$key}.".\n";
	}
	print LOGFILE "\n\n";
}

#checking that the input file contains only one sequence. To do this, I just calculate the number of lines starting with ">"
sub check_that_the_input_FASTA_file_contains_only_one_sequence()
{
	open STARTER, "< $starter_file_path";
	$number_of_lines_starting_with_greaterthan_sign=0;
	while(<STARTER>)
	{
		if($_=~/^>/)
		{
			$number_of_lines_starting_with_greaterthan_sign++;
		}
	}
	close(STARTER);
	
	if($number_of_lines_starting_with_greaterthan_sign!=1)
	{
		print LOGFILE "\nError. Only 1 sequence should be in the input FASTA file with the starter, while the file \"$starter_file_path\" has $number_of_lines_starting_with_greaterthan_sign.\n";
		exit();
	}
} 