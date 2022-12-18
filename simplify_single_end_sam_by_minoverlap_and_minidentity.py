#!/usr/bin/env python3
# coding=utf-8

"""
ВАЖНО: скрипт работает с sam-файлами, которые выдаются Minimap2. Будет ли он работать с sam-файлами, выдаваемыми другими программами, я не проверял.
ВАЖНО: скрипт смотрит только на первичное выравнивание (primary alignment). Вторичные он всегда выкидывает, даже если они проходят критерии. "Вторичные" это не когда был spliced-mapping, а именно вторичные для одного сегмента рида.

Этот скрипт берёт на вход sam-файл и оставляет в нём только те риды, которые выровнялись к референсу минимум на столько-то нуклеотидов с минимальным сходством по этому фрагменту таким-то. Например, оставляет только риды, от которых выровнялось минимум 5 килобаз со сходством по выровнявшейся части минимум 95%.

Скрипт сделан на основе скрипта /mnt/lustre/shelkmike/Work/Mine/Chloroplasts13/Fagopyrum/Esculentum/Cultivated_buckwheat/Nuclear/Calc/Phylogenetic_distance_by_genome_skimming/2_-_Extract_nuclear_reads/get_nuclear_reads.py . Если нужно будет сделать аналогичный скрипт для парноконцевых ридов, можно будет его делать на основе того же скрипта, там и для парноконцевых ридов есть часть.

Способ запуска:
simplify_single_end_sam_by_minoverlap_and_minidentity.py input_sam minoverlap minidentity output_sam

Пример:
simplify_single_end_sam_by_minoverlap_and_minidentity.py infile.sam 5000 0.95 outfile.sam
"""

import sys
import os
import subprocess
import re
import random
import datetime


f_infile = open(sys.argv[1], "r")
n_minimum_aligned_part_length = int(sys.argv[2])
n_minimum_sequence_similarity_in_the_aligned_part = float(sys.argv[3])
f_outfile = open(sys.argv[4], "w")
	
l_infile = list(f_infile)

for n_line_number in range(0, len(l_infile)):
	"""
	@SQ     SN:Dasha_plastid_genome__MT364821       LN:247391
	@SQ     SN:MT318703.1   LN:71837
	@SQ     SN:MT318705.1   LN:52654
	@SQ     SN:MT318701.1   LN:35904
	@PG     ID:minimap2     PN:minimap2     VN:2.22-r1101   CL:minimap2 -ax sr --sam-hit-only -o Dozhdik_1_mapping.sam -t 22 ./combined_reference.fasta /mnt/lustre/shelkmike/Work/Mine/Chloropla
	DRR046985.323736715     83      NODE_5755_length_64705_cov_54.22        25123   1       90M4S   =       25075   -138    TAAGATAGATACCGAGAAGTGTTTATATAAGAGCTTAGGAAGAAATGGTGGTATTCCAATGCAAATTCAATTTAATCGTGATAACTCAGTTAGG       EFDDEEEDDDDDDFFFFFHHHGHHIIJJJJJJIJJJJJJJJJJJJJJJJJJJJJIJJJJIJJJJJJJJJJJJJJJJJJJJJJIJJJJJHHHHHF  NM:i:2  ms:i:160        AS:i:160        nn:i:0  tp:A:P  cm:i:3       s1:i:54 s2:i:0  de:f:0.0222     rl:i:0
	DRR046985.323736715     163     NODE_5755_length_64705_cov_54.22        25075   1       101M    =       25123   138     TAATATTGAAGGACTAAAGGTTGTTATACTGATCCGATGTCACGTTAGTAAGATAGATACCGAGAAGTGTTTATATAAGAGCTTAGGAAGAAATGGTGGTA        CCCFFFFFHHHHHJJJJJJJIIJIJJJJJJJJJJJJJJIIJJJJHJJJGIFHJJJJJJJJJJIJIJI?EEHHHFFFFFFFEEEEEEDDDDDDDDDDACD;A   NM:i:2  ms:i:182        AS:i:182        nn:i:0       tp:A:P  cm:i:6  s1:i:54 s2:i:71 de:f:0.0198     rl:i:0
	DRR046985.323736729     83      NODE_450_length_581562_cov_89.22        392402  60      101M    =       392316  -187    CCCACCAGCTGGCATCCAGTAAATATCAATAATCAATAATAAAAAGAAGAAAGATAAAAACCAAGGAAATAAAGCAAAGCATTAGGGTTTGAGAGAGAGAG        DDDDEEEDEFFFFFGHHHHHJIIHHGIJJIJJIJJJJJJJJJJJJJJJJJIJJIHJJIIGJJJJJJJJIGJJJJJJJJJIJJJIJJJJHHHHHFFFFFCCC   NM:i:1  ms:i:196        AS:i:196        nn:i:0       tp:A:P  cm:i:15 s1:i:144        s2:i:0  de:f:0.0099     rl:i:0
	DRR046985.323736729     163     NODE_450_length_581562_cov_89.22        392316  56      101M    =       392402  187     TAAGCGAGGATGAATCCGTAACGAGGATGAGTAGGTACTCCCCGCTCCAAATCTGCCCTATTAACACCTCCTTGACTATGATGGTCCCCACCAGCTGGCAT        CCCFFFFFHHHHHJJJJJGHIIIFHIHIIGIBFHGBGIGIJJIJJJJHHHHEHFFFFFEEEEEEEDDDDDDDDDDDDDDCCEDCACDBDDDDDDDDCDDDC   NM:i:6  ms:i:150        AS:i:150        nn:i:0       tp:A:P  cm:i:3  s1:i:144        s2:i:0  de:f:0.0594     rl:i:0
	DRR046985.323736731     83      NODE_4496_length_64334_cov_61.55        62285   60      101M    =       62206   -180    ATACTCGGCTGTGAGTTTTTGTATGCACTTTGATATATGATATTTCAGTTTTATCGATGCTGGTTTTATATTTATTTGTGCCTTGTTGCAAGATTCGGTGA        CA8DBBDCDDDDDDDDDEEDFFFFCDFHHHHHHJJJIIJJJJJJJJJJJJJJJJJJJJJJJJJJJJIIJJIJJJJJJJJJJJJJJJJJHHHHHFFFFFCCC   NM:i:6  ms:i:142        AS:i:142        nn:i:0       tp:A:P  cm:i:5  s1:i:75 s2:i:44 de:f:0.0594     rl:i:0
	DRR046985.323736731     163     NODE_4496_length_64334_cov_61.55        62206   60      101M    =       62285   180     AGAGATCAGCCGAGAGTTTCCGTACATGTAGGATAGGTTACGGGGTTTGTATGTACAAATACATGGTTGTGTAGTTTTGATACTCGGCTGTGAGTTTTTGT        CCCFFFFFHHHHGJIJHIJJJIHIIJJJGGIJHGIIJHGIJJJJJBHIJJDGIJJJHHHHHHFFFFDFDDDCDDDDEDDDDEDDDDDDDDBDED@DDDDD@   NM:i:4  ms:i:166        AS:i:166        nn:i:0       tp:A:P  cm:i:4  s1:i:75 s2:i:44 de:f:0.0396     rl:i:0
	DRR046985.323736746     99      NODE_9488_length_233629_cov_107.58      200183  1       59S42M  =       200190  95      AGGATCCTTGCGTTCTCATGGCGCTCATCCTTGTATTCATTGAAGCGGGTCTTAAGAGCATCCATCCGCTACAAAAGTCCGGCTACCCATGGCGGCTGCTC        CCCFFFFFHHGHHJJJJIIJJJJIIJJJJJJJJHIJJJJJJJJJJJJIIBGGHHHHHHFFFFFEEEEDDDDDDDDADCDDDDDDDDDDDDDDDDDDDDDCD   NM:i:2  ms:i:64 AS:i:64 nn:i:0  tp:A:P  cm:i:1       s1:i:32 s2:i:0  de:f:0.0476     rl:i:0
	DRR046985.323736746     147     NODE_9488_length_233629_cov_107.58      200190  1       88M13S  =       200183  -95     CGCTACAAAAGTCCGGCTACCCATGGCGGCTGCTCATCCATCAATGCCGGTCGCTGCACCATCCCTTCCGGCGGTGCTCCCATCCAAGTTGGAGGAGCACC        DDDEDDDDDDB@DDDDDC?DDDDDDDDDDDDCDDDDDDDEDDDDDDBDDDFFGHHIG@JIIHJIIIHJJJJIJIIJJIIJJIGJJJJJHHHHHFFFFDBCB   NM:i:9  ms:i:90 AS:i:90 nn:i:0  tp:A:P  cm:i:2       s1:i:32 s2:i:0  de:f:0.1023     rl:i:0
	DRR046985.323736780     99      NODE_4822_length_24675_cov_84.67        7570    0       101M    =       7610    141     GTCCTAAGACTTCTAAACCCTCAGATCTCGCAAGATCATAACTCACAACGCAACTAGGATAAGTCAAACCTAACTCCAAAATCCACAAGATCATACAACCC        CCCFFFFFHHHHHJJJJJJJJJJJJIJJJIJJJJJJJJJJJJIJJJJJJJJJJJJJJIJJJGHHHHHFFFFFEEEEEECDDDDDDDDDDDDDDCDDDDDDD   NM:i:2  ms:i:182        AS:i:182        nn:i:0       tp:A:P  cm:i:7  s1:i:84 s2:i:0  de:f:0.0198     rl:i:0
	DRR046985.323736780     147     NODE_4822_length_24675_cov_84.67        7610    0       101M    =       7570    -141    ACTCACAACGCAACTAGGATAAGTCAAACCTAACTCCAAAATCCACAAGATCATACAACCCCCCCTCGAGCCTCAAAACGCTCACAAGATCAACCATATAA        A@:825@<>>>C@@>CCCCADCC@?CA:CA>A:C:DCDDA@:B@@DDDCECEEC@@70DBF<(HD6GIJJJIJJJIHGFIJJIIJJJIGHFHFFFFFFCCC   NM:i:2  ms:i:182        AS:i:182        nn:i:0       tp:A:P  cm:i:3  s1:i:84 s2:i:0  de:f:0.0198     rl:i:0
	"""
	
	#if n_line_number % 100000 == 0:
	#	f_logs.write("Analysed " + str(n_line_number) + " lines in the sam-file\n")
	
	#если в этой строке что-то есть, и она не начинается с "@", значит это строка про рид.
	if re.search(r"^[^\@]", l_infile[n_line_number]):
		o_regular_expression_results = re.search(r"^([^\t]+)\t([^\t]+)\t([^\t]+)\t[^\t]+\t[^\t]+\t([^\t]+)\t[^\t]+\t[^\t]+\t[^\t]+\t([^\t]+)\t([^\t]+).+NM:i:(\d+)", l_infile[n_line_number])
		if o_regular_expression_results:
			s_read_title = o_regular_expression_results.group(1)
			n_FLAG = int(o_regular_expression_results.group(2))
			s_contig_title = o_regular_expression_results.group(3)
			s_CIGAR = o_regular_expression_results.group(4)
			s_read_sequence = o_regular_expression_results.group(5)
			s_read_quality_string = o_regular_expression_results.group(6)
			n_edit_distance = int(o_regular_expression_results.group(7))
			
			#смотрю, установлен ли во FLAG бит 256. Если установлен, то это вторичное выравнивание, и я его не рассматриваю.
			s_binary_FLAG = str(bin(n_FLAG)) #"256" таким методом сконвертируется в 0b100000000
			if (n_FLAG < 256) or (int(s_binary_FLAG[-9:-8]) != 1): #смотрю, не является ли девятая цифра справа единицей. Кроме этого, вообще проверяю, что FLAG равен меньше 256, потому что в таком случае строка 0b100000000 будет короче, чем нужно для анализа "s_binary_FLAG[-9:-8]" — какие последствия будут в таком случае при расчёте выражения "int(s_binary_FLAG[-9:-8]) != 1" я не разбирался, поэтому для безопасности и использую дополнительно "n_FLAG < 256".
				#посчитаю суммарную длину выровнявшихся нуклеотидов. Для этого из длины рида вычту то, сколько было на краю рида нуклеотидов с soft clipping. Нуклеотиды c hard clipping и так не включены в последовательность рида.
				n_number_of_soft_clipped_bases = 0
				while re.search(r"(\d+)[S]", s_CIGAR):
					o_regular_expression_results = re.search(r"(\d+)[S]", s_CIGAR)
					n_number_of_soft_clipped_bases += int(o_regular_expression_results.group(1))
					s_CIGAR = re.sub(o_regular_expression_results.group(1) + r"[S]", "", s_CIGAR, 1)
				
				n_aligned_part_length = len(s_read_sequence) - n_number_of_soft_clipped_bases #сколько нуклеотидов от длины рида выровнялось
				
				#if n_aligned_part_length == 0:
					#print("Strange, n_aligned_part_length is 0 for " + l_infile[n_line_number])
				
				n_sequence_similarity_in_the_aligned_part = 1 - n_edit_distance / n_aligned_part_length
				
				#print("In the read " + s_read_title + " the alignment length is " + str(n_aligned_part_length) + ". The similarity between the aligned read part and the reference is " + str(n_sequence_similarity_in_the_aligned_part))
				
				if (n_aligned_part_length >= n_minimum_aligned_part_length) and (n_sequence_similarity_in_the_aligned_part >= n_minimum_sequence_similarity_in_the_aligned_part):
					f_outfile.write(l_infile[n_line_number])
	else:
		f_outfile.write(l_infile[n_line_number])