![Logo of Elloreas](https://gcf.fbb.msu.ru/shelkmike/Elloreas_logo/elloreas_logo.jpeg)

Elloreas is a targeted genome assembler, which means it takes a starter and then iteratively extends it with sequencing reads. This starter may be a contig from another assembly or just a random read. The principle of operation of Elloreas is similar to the principles of many other similar tools, like [TASR](https://github.com/bcgsc/TASR), [Mapsemsler2](http://colibread.inria.fr/software/mapsembler2/), [NOVOPlasty](https://github.com/ndierckx/NOVOPlasty) or [DistributedNucleatingAssembler](https://github.com/JGI-Bioinformatics/Kmernator). However, while those assemblers are designed for short reads, Elloreas was created for long reads, produced by sequencing machines of PacBio and Oxford Nanopore.

I created Elloreas to assemble the mitochondrial genome of Fagopyrum esculentum. A feature of this genome is that it consists of 10 circular chromosomes. Moreover, they can recombine, merging with each other or splitting into several smaller chromosomes. Multiple alternative forms of chromosomes coexist within a single plant, this is called "structural heteroplasmy". You can read more about this in the paper ["Mitochondrial genome of Fagopyrum esculentum and the genetic diversity of extranuclear genomes in buckwheat"](https://pubmed.ncbi.nlm.nih.gov/32408719/) where Elloreas is described. I thought that Elloreas may be useful for other bioinformaticians, this is why I deposited it on GitHub.
  
### How to run Elloreas.
A simple version, which will often be enough:  
`perl elloreas.pl --reads reads.fastq --starter starter.fasta`

A slightly more complicated version:  
`perl elloreas.pl --reads reads.fastq --starter starter.fasta --sequencing_technology oxford_nanopore --minimum_length_of_mapped_read_part 8000`

For a full list of parameters with suggestions on how to use them run  
`perl elloreas.pl --help`
  
### How Elloreas works.
A more detailed description may be found in manual.pdf, but in short it works like this:	
<img src="https://gcf.fbb.msu.ru/shelkmike/Elloreas_logo/scheme.jpeg?" width="80%" height="80%"/>
<br />

At each iteration Elloreas produces a list of alternative extensions.
  
<br />

### What is Elloreas used for?
Usually de novo genome assemblers cannot assemble each chromosome in a whole contig. They produce many contigs which are fragmented in places which are created by "forks" in the assembly graph. Forks are places with two or more alternative extensions, which usually arise because of repeats in genomes, like transposons. You can read more about forks, for example, [here](https://academic.oup.com/bib/article/10/4/354/299108).
Elloreas is useful to:
1) Understand which alternative extensions exist for your contig. Then you can choose the extension you like and assemble the genome by Elloreas further using this particular extension.
2) Assemble whole genomes. To do this, use a random read as a starter. Elloreas in not very fast - expect a rate of starter elongation of approximately 1000 bp/minute. Therefore, I don't recommend to use it for genomes larger than bacterial.
  
### Requirements:
Perl  
R  
minimap 2.19 or newer  
samtools  
BLAST+  
bedtools  
USEARCH  
MAFFT  
EMBOSS  
LASTZ  

I tested Elloreas with the following versions, but it should work with others as well: Perl 5.28, R 3.5.1, samtools 1.9, BLAST+ 2.9.0, USEARCH 11.0.667, MAFFT 7.402, EMBOSS 6.6.0, LASTZ 1.04. Paths to binaries of all these programs should be available through the environment variable $PATH.
  
### Output of Elloreas
The main files produced by Elloreas are:
1) elloreas_logs.txt - this file contains the basic information of how Elloreas worked. The first thing you may want to look at when Elloreas finishes is the end of this file - it provides the most important information of the Elloreas run.
2) history_of_alternative_extensions.txt - this file contains the list of alternative extensions present in each iteration of elongation. Only the top 5 most represented by reads extensions are listed there.
3) contig_from_the_final_iteration.fasta - this is the final contig produced by Elloreas.
4) coverage_plot_for_the_contig_from_the_final_iteration.jpeg and dotplot_for_the_contig_from_the_final_iteration.jpeg - these two figures represent a read coverage distribution and a dotplot for the final contig produced by Elloreas.
5) Files named like Iteration15/contig.it15.fasta - they contain the contig as it was after a particular iteration of elongation (after the 15th iteration in this example). Such contig files are formed at the end of each iteration.

You can read about less important files produced by Elloreas in manual.pdf.
   
<br />

### Caveats:
1) Long reads usually contain a large percent of sequencing errors, namely short indels and substitutions of one nitrogenous base by another. Some of these errors may persist in your assembly. Therefore, after you assemble the sequence, I highly recommend to polish it with programs such as [Pilon](https://github.com/broadinstitute/pilon/wiki) or [Racon](https://github.com/isovic/racon). 
Don't forget to include all sequences from your genome during polishing. Otherwise, there may be mistakes. For example, if you polish a plant's mitochondrial genome but don't add the plastid genome during polishing, some of plastid reads will map to the mitochondrial genome (because these two genomes have homologous sequences), which may produce a chimeric sequence. During the assembly stage, the plastid genome won't be a big problem, because those homologous sequences create forks, which are reported by Elloreas, so you can choose the mitochondrial extension and ignore the plastid one.
  
<br />

### Questions and answers:
1) How to cite Elloreas?<br />

Cite the paper ["Mitochondrial genome of Fagopyrum esculentum and the genetic diversity of extranuclear genomes in buckwheat"](https://pubmed.ncbi.nlm.nih.gov/32408719/) where Elloreas was first described.

2) What's the difference between this readme and the full version of the manual?<br />

The full version of the manual (manual.pdf) provided with each version of Elloreas is identical to the document you now read, except that the manual contains a more detailed description of the output files and the algorithm of Elloreas.

3) Elloreas extends the contig only in the 3' direction ("to the right"). Why not make it to extend the contig to the left too?<br />

One of the most important features of Elloreas is that it provides a user with a list of forks, when they are present. Making a user to examine forks on 3' and 5' ends simultaneously may confuse the user. If you need to extend the contig leftward, just take the reverse complement version of the contig, so left becomes right. If you're assembling a circular chromosome, you won't even need to do this, because during the extension the right end will join the left end.

4) Can Elloreas be used with short reads?<br />

Yes, it can. Minimap2, the read mapper which Elloreas uses, is capable of mapping short reads too. However, Elloreas will use paired-end reads as single-end reads. This is why I recommend to use dedicated targeted assemblers for short reads, like NOVOPlasty or TASR.

5) Long reads of PacBio and Oxford Nanopore usually contain many sequencing errors. Why not integrate a polishing process like in [Racon](https://github.com/lbcb-sci/racon) into Elloreas, so it polishes the contig after the assembly finishes?<br />

The problem with contig polishing by read mapping is that if you take only a single contig for polishing, reads from paralogs will map onto it too, so the polished sequence will become a chimera. This is why it is important to add all contigs of your genome to the contig created by Elloreas and perform polishing only after this, thus reducing the possibility of false mappings.

6) There was a fork at some iteration. Elloreas chose one of alternative extensions, but I want to try another one. What should I do?<br />

In the file history_of_alternative_extensions.txt find the respective iteration. You will see a list of alternative extensions for this iteration. Elloreas always chooses the extension supported by the highest number of reads. Take the file contig_from_the_final_iteration.fasta, remove everything starting with this most-supported extension and instead paste the extension you want to try. Then, use this newly formed contig as a starter.

7) I'm using a read as a starter but Elloreas stops after the first iteration. Why is this so?<br />

The problem is that reads may be too dissimilar to each other. If sequencing errors make each of two reads 85% similar to the genome, then, if the sequencing errors were random, these two reads will be approximately 72% similar to each other (0.85\*0.85=0.7225). This is why reads hardly map to a starter which is also a read. To deal with this, I recommend to decrease the value of the "--minimum_read_similarity" parameter.

8) Where can I ask questions about Elloreas?<br />

You can use the "Issues" section on GitHub (https://github.com/shelkmike/Elloreas/issues)

9) How fast is Elloreas?<br />

I tested it using 22 CPU cores for several datasets each containing about a million reads and the rate of contig elongation was on the order of 1 kilobase per minute.

10) Will Elloreas run on Windows or Mac?<br />

I didn't test it on Mac, but I think it will work if you can install all the required programs. To run Elloreas on Windows use a virtual Linux machine, for example Ubuntu run through VirtualBox (https://www.osboxes.org/ubuntu/).

11) Many de novo assemblers create assembly graphs. Why can't they be used to detect forks in the assembly instead of Elloreas?<br />

You can use them. For example, Spades produces FASTG files with assembly graphs which can be visualised with Bandage. However, in my experience, de novo assemblers often don't report all forks. Probably, this is because there are some inner requirements of how represented by reads a branch should be to report this branch in a FASTG file. Elloreas, instead, produces a more detailed list of alternative extensions arising during the assembly process.

12) Well, finally, will Elloreas be useful for me?<br />

I don't know. It definitely was useful for ME. You may give it a try.
