About:
------
The program gFetch.pl accepts a text file (containing desired taxa in each row) from user as "input" and download the closest NCBI RefSeq/GenBank sequences in the user given output folder.

Commandline parameters:
-----------------------
gFetch.pl -i input_file -o output_directory 

Optional Input:
---------------
	-use_refseq 0/1 		[Default is 1; The program will download RefSeq sequences]
	-use_genbank 0/1 		[Default is 0; Setting it 1 will instruct the program to download GenBank sequences]

	-T number_of_threads		[Default is 4;] 

	-taxa_mapping_table table.txt 	[Unless a file is provided the program will create a file /tmp/query_and_target_taxa_with_downloaded_sequence.txt which will show the associated Taxa of downloaded sequences]

	-input_is_accession_list 1	[This option must be provided if the input file contains a list of NCBI nucleotide accessions]


Dependencies:
-------------
Fetching sequences using NCBI accession number requires the Entrez Programming Utilities (E-utilities) tools pre-installed in the system. Please refer to https://www.ncbi.nlm.nih.gov/books/NBK179288/
for detail instruction of how to install the tools.

For Machintosh users: the "wget" tool is needed to be installed.


Note:
-----
Please note that the input text will be used in a "case-sensitive" manner.

This tool is developed and tested in Fedora 22 server and should run in all *nix operating systems.
