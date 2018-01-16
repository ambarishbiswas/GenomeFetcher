#!/usr/bin/perl -w
use strict;
use Cwd 'abs_path';

use FindBin qw($Bin);
use lib "$Bin/lib1";                # where was script installed?



#------ the user input parameters -----------

my $input_file="";
my $output_file="";
our $temp_dir="";

my $number_of_threads=8;
my $use_refseq=1;
my $use_genbank=0;

my $input_is_accession_list=0;

my $update_reference_tables=0;
my $clean_temp_files=0;


#-------Reference files ---------------------

our $table_with_RefSeq_taxa_and_associated_sequence_links="$Bin/REF_FILES/table_with_RefSeq_taxa_and_associated_sequence_links.txt";
our $table_with_GenBank_taxa_and_associated_sequence_links="$Bin/REF_FILES/table_with_GenBank_taxa_and_associated_sequence_links.txt";





#-----------------------------------------------------
for(my $i=0;$i<=$#ARGV;$i++)
	{
						
		#--------------------------- input options ------------------------------------------------------------
		if($ARGV[$i]=~/-i$/)
			{
				$input_file=$ARGV[$i+1];					
			}
		elsif($ARGV[$i]=~/-input_is_accession_list$/i)
			{
				$input_is_accession_list=$ARGV[$i+1];					
			}
		elsif($ARGV[$i]=~/-o$/)
			{
				$output_file=$ARGV[$i+1];					
			}	
		elsif($ARGV[$i]=~/-temp_dir$/)
			{
				$temp_dir=$ARGV[$i+1];					
			}											
		elsif($ARGV[$i]=~/-T$/)
			{
				$number_of_threads=$ARGV[$i+1];					
			}
		elsif($ARGV[$i]=~/-use_refseq$/i)
			{
				$use_refseq=$ARGV[$i+1];					
			}
		elsif($ARGV[$i]=~/-use_genbank$/i)
			{
				$use_genbank=$ARGV[$i+1];					
			}
		elsif($ARGV[$i]=~/-update_reference_tables$/i)
			{
				$update_reference_tables=$ARGV[$i+1];					
			}
		elsif($ARGV[$i]=~/-clean_temp_files$/i)
			{
				$clean_temp_files=$ARGV[$i+1];					
			}		
	}
#--- create temp directory unless provided
if($update_reference_tables==1 or $input_file=~/\S/ or $output_file=~/\S/ or $temp_dir=~/\S/)
	{
		unless(-d $temp_dir)
			{
				if($temp_dir!~/\S/)
					{
						$temp_dir="/tmp/".&get_unique_id();
					}
					
				system("mkdir -p $temp_dir");
				
				print "\nAll intermediate files and sequence(s) will be downloaded in this folder: \n\t$temp_dir\n";
			}
	}	
	
if($update_reference_tables==1)
	{
		&check_dependencies($update_reference_tables);
		$update_reference_tables=0;
		print "\nAll reference files are updated.\n\n"; 
		
		if($input_file!~/\S/ and $output_file!~/\S/)
			{
				exit;
			}	
	}

if($input_file!~/\S/ or $output_file!~/\S/)
	{
		print "Error: please provide the input and out put files in this manner:\n\t$0 -i input_file -o output_directory [Optional: -input_is_accession_list 0/1 -use_refseq 1/0 -use_genbank 1/0 -T number_of_threads -update_reference_tables 0/1 ]\n\n"; exit;
	}
	
	
#---- do not set the number of threads too high; NCBI terminates connection--	
if($number_of_threads > 8){$number_of_threads=8;}
use Parallel::ForkManager;
my $pm = new Parallel::ForkManager($number_of_threads);


#--- set the other two files --------------------------------
my $accession_and_taxID_table	=$output_file.".accession_and_taxID_table.txt";
my $taxa_mapping_table			=$output_file.".query_term_and_matched_Taxa.txt";
#-- ---------------------------------------------------


if($input_is_accession_list==1)
	{
		&get_sequence_using_accessions($pm,$input_file);
	}
else{	
	
		#---- first check if the reference files are present--
		&check_dependencies($update_reference_tables);
		
		#---- load the taxa and ftp for the selected group
		my %hash_of_ftp_and_associated_taxID;
		my %hash_of_full_taxa_and_ftps;
		my %hash_of_lvl0_taxa_full_taxa_and_ftps;
		my %hash_of_lvl1_taxa_full_taxa_and_ftps;
		my %hash_of_lvl2_taxa_full_taxa_and_ftps;

		if($use_refseq==1)
			{
				print "Loading RefSeq taxa and ftp files..\n";
				my @arr_rows=`cat $table_with_RefSeq_taxa_and_associated_sequence_links >&1`;
				foreach my $row(@arr_rows)
					{
						chomp $row; $row=~s/\r//;
						
						my($taxa,$ftps,$taxID)=split('\t',$row);
						
						my @arr_t1=split(';',$taxa);
						my @arr_taxa;
						my $count=0;
						foreach my $name(@arr_t1)
							{
								if($name=~/\s/)
									{
										#-- skip
									}
								else{
										push(@arr_taxa,$name);
										$count++;
									}	
								if($count>3){last;}	
							}
						#my @arr_taxa=split(';',$taxa);
						
						
						$hash_of_full_taxa_and_ftps{$taxa}=$ftps;
						
						if($arr_taxa[0] and $arr_taxa[0]=~/\S/)
							{
								$hash_of_lvl0_taxa_full_taxa_and_ftps{$arr_taxa[0]}{$taxa}=$ftps;
							}
						if($arr_taxa[1] and $arr_taxa[1]=~/\S/)
							{
								$hash_of_lvl1_taxa_full_taxa_and_ftps{$arr_taxa[0]}{$arr_taxa[1]}{$taxa}=$ftps;
							}
						if($arr_taxa[2] and $arr_taxa[2]=~/\S/)
							{
								$hash_of_lvl2_taxa_full_taxa_and_ftps{$arr_taxa[0]}{$arr_taxa[1]}{$arr_taxa[2]}{$taxa}=$ftps;
							}
						
						
						my @arr_ftp=split(';',$ftps);
						foreach my $ftp(@arr_ftp)
							{
								$hash_of_ftp_and_associated_taxID{$ftp}=$taxID;
							}					
							
					}
				
			}

		if($use_genbank==1)
			{
				print "Loading GenBank taxa and ftp files..\n";
				my @arr_rows=`cat $table_with_GenBank_taxa_and_associated_sequence_links >&1`;
				foreach my $row(@arr_rows)
					{
						chomp $row; $row=~s/\r//;
						
						my($taxa,$ftps,$taxID)=split('\t',$row);
						
						my @arr_t1=split(';',$taxa);
						my @arr_taxa;
						my $count=0;
						foreach my $name(@arr_t1)
							{
								if($name=~/\s/)
									{
										#-- skip
									}
								else{
										push(@arr_taxa,$name);
										$count++;
									}	
								if($count>3){last;}	
							}
						#my @arr_taxa=split(';',$taxa);
						
						$hash_of_full_taxa_and_ftps{$taxa}=$ftps;
						
						if($arr_taxa[0] and $arr_taxa[0]=~/\S/)
							{
								$hash_of_lvl0_taxa_full_taxa_and_ftps{$arr_taxa[0]}{$taxa}=$ftps;
							}
						if($arr_taxa[1] and $arr_taxa[1]=~/\S/)
							{
								$hash_of_lvl1_taxa_full_taxa_and_ftps{$arr_taxa[0]}{$arr_taxa[1]}{$taxa}=$ftps;
							}
						if($arr_taxa[2] and $arr_taxa[2]=~/\S/)
							{
								$hash_of_lvl2_taxa_full_taxa_and_ftps{$arr_taxa[0]}{$arr_taxa[1]}{$arr_taxa[2]}{$taxa}=$ftps;
							}
						
						
						
						
						my @arr_ftp=split(';',$ftps);
						foreach my $ftp(@arr_ftp)
							{
								$hash_of_ftp_and_associated_taxID{$ftp}=$taxID;
							}					
							
					}
				
			}



			
		#--- process the input file ---
		print "\n\tProcessing user inputted file..\n";
		my @arr_inputted_taxa=`grep -v '^#' $input_file | cut -f1 >&1`;

		my %hash_of_inputted_row_and_cleaned_taxa;
		my %hash_of_inputted_taxa;

		my @arr_inputted_original_rows;
		foreach my $row(@arr_inputted_taxa)
			{
				chomp $row; $row=~s/\r//; 
				
				my $original_row=$row;		
				
				if($row=~/\w{1}__/)
					{
						$row=~s/\w{1}__//g;
					}
				
				#--- remove all underscores with space
				$row=~s/_/ /g;
				if($hash_of_inputted_taxa{$row})
					{
						$hash_of_inputted_taxa{$row}=$hash_of_inputted_taxa{$row}+1;
					}
				else{
						$hash_of_inputted_taxa{$row}=1;
					}	
					
				$hash_of_inputted_row_and_cleaned_taxa{$original_row}=$row;		
				push(@arr_inputted_original_rows,$original_row);	
				
					
			}

		my %hash_of_qry_taxa_and_ftps;
		my %hash_of_qry_taxa_and_nof_required_ftps;

		foreach my $row(sort{length($b)<=>length($a)} keys %hash_of_inputted_taxa)
			{
				chomp $row; $row=~s/\r//;#$row=$row.";";
				
				
				#print "Processing: $row\n";
				
				if($row=~/\w{1}__/)
					{
						$row=~s/\w{1}__//g;
					}
				
				#--- remove all underscores with space
				$row=~s/_/ /g;
				
				#----
				my @arr_r1;
				if($row=~/;/)
					{		
						@arr_r1=split(';',$row);
					}
				else{
						push(@arr_r1,$row);
					}		
					
				my @arr_taxa;
				foreach my $name(@arr_r1)
					{
						if($name and $name=~/\w+/)
							{
								push(@arr_taxa,$name);
							}
					}
				
				my @arr_ftps;
				my %hash_of_ftps;
				my $match_found=0;
				
				
				#--- check if only one name provided and its Root
				if($row eq 'Root' or $row eq 'root')
					{
						#-- check how many sequences 
						my $q_taxa=$hash_of_inputted_row_and_cleaned_taxa{$row};		
						my $nof_required_ftps=$hash_of_inputted_taxa{$q_taxa};
						
						my $ftp_count=0;
						foreach my $taxa(keys %hash_of_full_taxa_and_ftps)
							{
								my $ftp=$hash_of_full_taxa_and_ftps{$taxa};
								$hash_of_qry_taxa_and_ftps{$row}{$ftp}=$taxa;
								$match_found=1;
								$ftp_count++;
								if($ftp_count>=$nof_required_ftps){last;}
							}
							
					}		
				
				#-- check the lvl2 hash
				if($match_found==3 and $arr_taxa[2] and $arr_taxa[2]=~/\S/ and $hash_of_lvl2_taxa_full_taxa_and_ftps{$arr_taxa[0]}{$arr_taxa[1]}{$arr_taxa[2]})
					{	
						#print "\t2:$row\n";	
						
						my %hash_of_target_taxas;	
						my %hash_of_matched_taxas;
							
						foreach my $taxa(keys %{$hash_of_lvl2_taxa_full_taxa_and_ftps{$arr_taxa[0]}{$arr_taxa[1]}{$arr_taxa[2]}})
							{						
								#print "\t$taxa\n";
								my $ftps=$hash_of_lvl2_taxa_full_taxa_and_ftps{$arr_taxa[0]}{$arr_taxa[1]}{$arr_taxa[2]}{$taxa};
								$hash_of_target_taxas{$taxa}=$ftps;				
								
							}
						
						my($highest_score)=&identify_closest_taxa_and_get_ftps($row,\%hash_of_target_taxas,\%hash_of_matched_taxas,\%hash_of_qry_taxa_and_ftps);	
						
						if($highest_score==0)
							{
								#print "\t2: Error: no match found for $row\n";next;
							}
						else{
								$match_found++;	
									
							}
						
						#foreach my $matched_taxa(keys %hash_of_matched_taxas)
						#	{
						#		#print "\tHS:$highest_score: $matched_taxa\t$hash_of_matched_taxas{$matched_taxa}\n";
						#	}		
						#print "\n";						
					}
				
				if($match_found==0 and $arr_taxa[1] and $arr_taxa[1]=~/\S/ and $hash_of_lvl1_taxa_full_taxa_and_ftps{$arr_taxa[0]}{$arr_taxa[1]})
					{		
						#print "\t1:$row\n";
						
						my %hash_of_target_taxas;	
						my %hash_of_matched_taxas;		
						foreach my $taxa(keys %{$hash_of_lvl1_taxa_full_taxa_and_ftps{$arr_taxa[0]}{$arr_taxa[1]}})
							{
								#print "\t$taxa\n";
								my $ftps=$hash_of_lvl1_taxa_full_taxa_and_ftps{$arr_taxa[0]}{$arr_taxa[1]}{$taxa};
								$hash_of_target_taxas{$taxa}=$ftps;	
							}
						
						my($highest_score)=&identify_closest_taxa_and_get_ftps($row,\%hash_of_target_taxas,\%hash_of_matched_taxas,\%hash_of_qry_taxa_and_ftps);	
						
						if($highest_score==0)
							{
								#print "\tError: no match found for $row\n";next;
							}
						else{
								$match_found++;		
							}
							
						#foreach my $matched_taxa(keys %hash_of_matched_taxas)
						#	{
						#		print "\tHS:$highest_score: $matched_taxa\t$hash_of_matched_taxas{$matched_taxa}\n";
						#	}		
						#print "\n";								
					}	
				
				if($match_found==0 and $arr_taxa[0] and $arr_taxa[0]=~/\S/ and $hash_of_lvl0_taxa_full_taxa_and_ftps{$arr_taxa[0]})
					{		
						#print "\t0:$row\n";
						
						my %hash_of_target_taxas;	
						my %hash_of_matched_taxas;		
						foreach my $taxa(keys %{$hash_of_lvl0_taxa_full_taxa_and_ftps{$arr_taxa[0]}})
							{						
								my $ftps=$hash_of_lvl0_taxa_full_taxa_and_ftps{$arr_taxa[0]}{$taxa};
								$hash_of_target_taxas{$taxa}=$ftps;	
							}	
						my($highest_score)=&identify_closest_taxa_and_get_ftps($row,\%hash_of_target_taxas,\%hash_of_matched_taxas,\%hash_of_qry_taxa_and_ftps);	
						
						if($highest_score==0)
							{
								#print "\tError: no match found for $row\n";#next;
							}
						else{
								$match_found++;		
							}
						#print "\tScore:$highest_score\n";
						#foreach my $matched_taxa(keys %hash_of_matched_taxas)
						#	{
						#		print "\t$matched_taxa\t$hash_of_matched_taxas{$matched_taxa}\n";
						#	}		
						#print "\n";						
					}
					
				if($match_found==0)
					{
						#--- search all the taxa for the closest match
						
						my %hash_of_matched_taxas;	
						my($highest_score)=&identify_closest_taxa_and_get_ftps($row,\%hash_of_full_taxa_and_ftps,\%hash_of_matched_taxas,\%hash_of_qry_taxa_and_ftps);	
						
						if($highest_score==0)
							{
								#print "\t-1:Error: no match found for $row\n";#next;
							}
						else{
								$match_found++;		
							}	
						
						my $total_matched_taxas=keys %hash_of_matched_taxas;
						
						#print "\tScore:$highest_score\t\$total_matched_taxas=$total_matched_taxas\n";
					

					}
					
					
				if($match_found==0)
					{
						print "Error: no match found for $row\n";
					}					
			}



		#---- now print the rows and a randomly selected ftplink
		my %hash_of_used_ftps;
		my $count=0;
		open(TAB,">$taxa_mapping_table");
		foreach my $original_row(@arr_inputted_original_rows)
			{
				$count++;
				
				my $q_taxa=$hash_of_inputted_row_and_cleaned_taxa{$original_row};		
				my $nof_required_ftps=$hash_of_inputted_taxa{$q_taxa};	
				
				my $ftp_count=0;
				my $matched_taxa="NA";
				my $matched_ftp="NA";
				
				
				#print "$count: $original_row\n";
				
				if($hash_of_qry_taxa_and_ftps{$q_taxa})
					{
						foreach my $ftp(keys %{$hash_of_qry_taxa_and_ftps{$q_taxa}})
							{
								my $t_taxa=$hash_of_qry_taxa_and_ftps{$q_taxa}{$ftp};
								
								if(not $hash_of_used_ftps{$ftp})
									{
										#print "\t$t_taxa\n";
										#print "\t$ftp\n";
										#print "\n";
										
										$matched_taxa=$t_taxa;
										$matched_ftp=$ftp;
										
										#print TAB "$original_row\t$t_taxa\t$ftp\n";
										
										$hash_of_used_ftps{$ftp}="$q_taxa\t$t_taxa\t$ftp";						
										$ftp_count++;
										
										last;
									}
								if($ftp_count >= $nof_required_ftps){last;}	
							}
					}	
				
				if($ftp_count==0)
					{		
						print "\tWarning: no sequence found for $original_row\n";
					}
				
				print TAB "$original_row\t$matched_taxa\t$matched_ftp\n";
			}
		close(TAB);

		#----------------------------------------------------


		my $total_files=keys %hash_of_used_ftps;
		my $index=0;

		foreach my $ftp(keys %hash_of_used_ftps)
			{
				chomp $ftp; $ftp=~s/\r//; $ftp=~s/;$//;
				$index++;
				my $remaining=$total_files-$index;
				
				print "\tDownloading $remaining/$total_files files..\n";
				
				my @arr_f=split('\/',$ftp);
				my $only_filename=$arr_f[$#arr_f];
				
				select(undef,undef,undef,0.5);
				
				$pm->start and next;
				system("wget -q '$ftp' -O '$temp_dir/$only_filename'");
				
				my $taxID=$hash_of_ftp_and_associated_taxID{$ftp};
				
				#--- extract the accession and make a table that will have accession and taxID
				my $ret=`zcat $temp_dir/$only_filename | grep '^>' | sed 's/>//' | awk -F '.' '{print \$1"\t"$taxID}' >$temp_dir/$only_filename.acc_and_taxID.txt`;
				#--
				$pm->finish;
			}
		$pm->wait_all_children;

}



#--- now zcat the downloaded sequences in the output file
my $ret2=`zcat $temp_dir/*.fna.gz >$output_file`; 

#--- also cat the accession and taxIDs
my $ret3=`cat $temp_dir/*.acc_and_taxID.txt >$accession_and_taxID_table`;

print "\nThe following files are created:\n";
print "\t1. $output_file                                [Contains the downloaded sequences]\n"; 
print "\t2. $accession_and_taxID_table  [Contains the RefSeq/GenBank sequence accession and associated taxon ID]\n"; 


if($input_is_accession_list==0)
	{
		print "\t3. $taxa_mapping_table [Table with query row/term, matched taxa and sequence_source_ftp_link]\n"; 
	}
#-- remove the tmp directory
if($clean_temp_files==1)
	{
		print "\n\tCleaning temporary files..\n";
		system("rm -rf .$temp_dir");
	}
else{	
		print "\nThe following directory contains all the intermediate files:\n";
		print "\t$temp_dir\n\n"; 
	}	


print "\nJob completed successfully.\n";

exit;













################## subs

sub check_dependencies()
	{
		my $update_reference_tables=shift(@_);
		

		#--- create the REF_FILES folder unless exist
		unless(-d "$Bin/REF_FILES")
			{
				system("mkdir -p $Bin/REF_FILES");
			}


		#--- now process the taxdump files and create the  exist table_with_RefSeq_taxa_and_associated_sequence_links and table_with_GenBank_taxa_and_associated_sequence_links tables
		if(not -e $table_with_RefSeq_taxa_and_associated_sequence_links or not -e $table_with_GenBank_taxa_and_associated_sequence_links or $update_reference_tables!=0)
			{			
				my %hash_of_taxID_and_name;
				my %hash_of_name_and_taxID;
				
				my %hash_of_taxID_and_parentTaxID;
				my %hash_of_parentTaxID_and_childTaxID;
				
				#----------------------------------------------------------------
				
					
				my $names_dmp_file="$temp_dir/names.dmp";
				my $nodes_dmp_file="$temp_dir/nodes.dmp";
																							#my $accession2taxid_file="nucl_gb.accession2taxid";
				my $assembly_summary_refseq="$temp_dir/assembly_summary_refseq.txt"; 	#ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
				my $assembly_summary_genbank="$temp_dir/assembly_summary_genbank.txt";

										
				#--- create a tmp directory and download the files in it
				unless(-d "$temp_dir"){system("mkdir -p $temp_dir");}	
											
				#--- wget the latest version of taxdumpfiles, RefSeq and Genbank assembly summary files --
				print "\nDownloading the following files ..\n";
				print "\tlatest NCBI taxonomy files ..\n";
				system("wget -q --timestamping ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -O $temp_dir/taxdump.tar.gz ");
				
				#--- wget the latest version of RefSeq assembly summary
				print "\tlatest RefSeq assembly summary file ..\n";
				system("wget -q --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt -O $temp_dir/assembly_summary_refseq.txt ");
				
				#--- wget the latest version of GenBank assembly summary
				print "\tlatest GenBank assembly summary file ..\n";
				system("wget -q --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt -O $temp_dir/assembly_summary_genbank.txt");
				
				
				
				
				
						
				#-- extract it
				if(-e "$temp_dir/taxdump.tar.gz")
					{ 
						#print "\tExtracting files ..\n";
						print "\nProcessing downloaded files ..\n";						
						print "\tNCBI taxonomy files ..\n";
						
						system("tar -C $temp_dir -zxvf $temp_dir/taxdump.tar.gz >/dev/null");
						
						#my $pwd=`pwd >&1`; 	print "***** $pwd \n\n";
						
						#-- get only the "scientific name" from names.dmp
						system("grep 'scientific name' $temp_dir/names.dmp >$temp_dir/names.dmp2");
						
						#--- copy the nodes.dmp file
						system("cp $temp_dir/names.dmp2 $names_dmp_file");
						
						#--- now process the files and load hashes
						
						open(RD,"$names_dmp_file");
						while(my $row=<RD>)
							{
								my @arr_r1=split('\t\|\t',$row);
								my $taxID=$arr_r1[0]; $taxID=int($taxID);
								my $name=$arr_r1[1];
												
								if($arr_r1[2] and $arr_r1[2]=~/</ and $arr_r1[2]=~/>/ ) #1987911	|	Zeugodacus	|	Zeugodacus <subgenus>	|	scientific name	|
									{
										#$name=$arr_r1[2];
									}									
								unless($name eq 'root')
									{
										$hash_of_taxID_and_name{$taxID}=$name;
										$hash_of_name_and_taxID{$name}{$taxID}=1;												
									}				
							}
						close(RD);
						
						#----------------------------------------------------------------
						#print "Loading $nodes_dmp_file\n";
						
						open(RD,"$nodes_dmp_file");
						while(my $row=<RD>)
							{
								my @arr_r1=split('\t\|\t',$row);
								
								my $taxID=$arr_r1[0]; $taxID=int($taxID);
								my $parent_taxID=$arr_r1[1]; $parent_taxID=int($parent_taxID);			
								
								if($parent_taxID >1)
									{
										$hash_of_taxID_and_parentTaxID{$taxID}=$parent_taxID;
										$hash_of_parentTaxID_and_childTaxID{$parent_taxID}=$taxID;
									}				
							}
						close(RD);
						
					}
				else{
						print "\n\nError: Couldn't download latest taxonomy files. Check network connection.\n\n"; exit;
					}						
				

				#-----------------------------------------------
				
				if(-e "$assembly_summary_refseq")
					{ 
						print "\tRefSeq assembly summary ...\n";
						my %hash_of_refseq_taxID_and_sequence_ftp_path;
						my %hash_of_refseq_ftp_and_taxIDs;
						
						my @arr_refseq_taxID_and_file_path=`awk -F "\t" '{if(\$20 ~ "ftp:"){print \$6"\t"\$20}}' $assembly_summary_refseq | awk -F "/" '{print \$0"/"\$NF"_genomic.fna.gz"}' >&1`;
						foreach my $row(@arr_refseq_taxID_and_file_path)
							{
								chomp $row; $row=~s/\r//;
								my @arr_r1=split('\t',$row);
								my $taxID=$arr_r1[0]; $taxID=int($taxID);
								my $ftp=$arr_r1[1];
								
								$hash_of_refseq_ftp_and_taxIDs{$ftp}=$taxID;
								$hash_of_refseq_taxID_and_sequence_ftp_path{$taxID}{$ftp}=1;				
							}
						
						#--- now create a table that will hold the name of each taxa level and all availabe child taxIDs with sequence
								
						#print "\tCreating $table_with_RefSeq_taxa_and_associated_sequence_links ..\n";	
						open(TAB2,">$table_with_RefSeq_taxa_and_associated_sequence_links");		
						foreach my $tax_name(sort keys %hash_of_name_and_taxID)
							{
								if($tax_name eq "cellular organisms"){next;}
								
								foreach my $taxID(keys %{$hash_of_name_and_taxID{$tax_name}})
									{
											
										#my $taxID=$hash_of_name_and_taxID{$tax_name};
										
										my %hash_of_current_taxIDs;				
										$hash_of_current_taxIDs{$taxID}=1;
										
										my %hash_of_parentNames;
										
										my @arr_taxa;
										unshift(@arr_taxa,$tax_name);
										
										
										#-------
										if($hash_of_taxID_and_parentTaxID{$taxID})
											{				
												#--- this block will get all the associated downstream sequences joined
												my $last_parentTaxID=$taxID;
												while($last_parentTaxID >1)
													{	
														if(not $hash_of_taxID_and_parentTaxID{$last_parentTaxID}){last;}
																											
														my $parent_taxID=$hash_of_taxID_and_parentTaxID{$last_parentTaxID};
														my $parent_name=$hash_of_taxID_and_name{$parent_taxID};		
														
														if($parent_name=~/\S+/ and $parent_name!~/cellular organisms/)
															{
																unshift(@arr_taxa,$parent_name);
																$last_parentTaxID=$parent_taxID;
															}
														else{
																$last_parentTaxID=0;
																last;
															}								
													}							
											}
										#-------

										my $row=join(";",@arr_taxa);
										$row=$row.";";
										
										#--- get the ftps
										my @arr_ftps;
										if($hash_of_refseq_taxID_and_sequence_ftp_path{$taxID})
											{
												foreach my $ftp(keys %{$hash_of_refseq_taxID_and_sequence_ftp_path{$taxID}})
													{
														push(@arr_ftps,$ftp);
													}						
											}
											
										my $associated_sequences="NA";
										if($arr_ftps[0])
											{
												$associated_sequences=join(";",@arr_ftps);
												$associated_sequences=$associated_sequences.";";
											}
										if($associated_sequences ne "NA")
											{		
												print TAB2 "$row\t$associated_sequences\t$taxID\n";
											}	
									}	
														
							} 
						close(TAB2);
					}
				else{
						print "\n\nError: Couldn't download latest RefSeq assembly summary file. Check network connection.\n\n"; exit;
					}
				#--------------

				
				
				
				if(-e "$assembly_summary_genbank")
					{ 
						print "\tGenbank assembly summary ...\n";
						my %hash_of_genbank_taxID_and_sequence_ftp_path;
						my %hash_of_genbank_ftp_and_taxIDs;
						
						my @arr_genbank_taxID_and_file_path=`awk -F "\t" '{if(\$20 ~ "ftp:"){print \$6"\t"\$20}}' $assembly_summary_genbank | awk -F "/" '{print \$0"/"\$NF"_genomic.fna.gz"}' >&1`;
						foreach my $row(@arr_genbank_taxID_and_file_path)
							{
								chomp $row; $row=~s/\r//;
								my @arr_r1=split('\t',$row);
								my $taxID=$arr_r1[0]; $taxID=int($taxID);
								my $ftp=$arr_r1[1];
								
								$hash_of_genbank_ftp_and_taxIDs{$ftp}=$taxID;
								$hash_of_genbank_taxID_and_sequence_ftp_path{$taxID}{$ftp}=1;				
							}
						
						#--- now create a table that will hold the name of each taxa level and all availabe child taxIDs with sequence
								
						#print "\tCreating $table_with_GenBank_taxa_and_associated_sequence_links ..\n";	
						open(TAB3,">$table_with_GenBank_taxa_and_associated_sequence_links");		
						foreach my $tax_name(sort keys %hash_of_name_and_taxID)
							{
								if($tax_name eq "cellular organisms"){next;}
								
								foreach my $taxID(keys %{$hash_of_name_and_taxID{$tax_name}})
									{
											
										#my $taxID=$hash_of_name_and_taxID{$tax_name};
										
										my %hash_of_current_taxIDs;				
										$hash_of_current_taxIDs{$taxID}=1;
										
										my %hash_of_parentNames;
										
										my @arr_taxa;
										unshift(@arr_taxa,$tax_name);
										
										
										#-------
										if($hash_of_taxID_and_parentTaxID{$taxID})
											{				
												#--- this block will get all the associated downstream sequences joined
												my $last_parentTaxID=$taxID;
												while($last_parentTaxID >1)
													{	
														if(not $hash_of_taxID_and_parentTaxID{$last_parentTaxID}){last;}
																											
														my $parent_taxID=$hash_of_taxID_and_parentTaxID{$last_parentTaxID};
														my $parent_name=$hash_of_taxID_and_name{$parent_taxID};		
														
														if($parent_name=~/\S+/ and $parent_name!~/cellular organisms/)
															{
																unshift(@arr_taxa,$parent_name);
																$last_parentTaxID=$parent_taxID;
															}
														else{
																$last_parentTaxID=0;
																last;
															}								
													}							
											}
										#-------

										my $row=join(";",@arr_taxa);
										$row=$row.";";
										
										#--- get the ftps
										my @arr_ftps;
										if($hash_of_genbank_taxID_and_sequence_ftp_path{$taxID})
											{
												foreach my $ftp(keys %{$hash_of_genbank_taxID_and_sequence_ftp_path{$taxID}})
													{
														push(@arr_ftps,$ftp);
													}						
											}
											
										my $associated_sequences="NA";
										if($arr_ftps[0])
											{
												$associated_sequences=join(";",@arr_ftps);
												$associated_sequences=$associated_sequences.";";
											}
										if($associated_sequences ne "NA")
											{		
												print TAB3 "$row\t$associated_sequences\t$taxID\n";
											}	
									}	
														
							} 
						close(TAB3);
				
					}
				
				else{
						print "\n\nError: Couldn't download latest RefSeq assembly summary file. Check network connection.\n\n"; exit;
					}
				#------ create a table with taxonID and full taxa
				my $table_with_taxID_and_full_taxa="$Bin/REF_FILES/table_with_taxID_and_full_taxa.txt";
				my $skip_construction_of_full_ncbi_taxa=1;
				if($skip_construction_of_full_ncbi_taxa==0)
					{
						#--- now write the full taxa for all available taxIDs----
						print "Creating $table_with_taxID_and_full_taxa\n";
						open(TAB,">$table_with_taxID_and_full_taxa");
						my %hash_of_already_created_taxaID_and_taxa;
						foreach my $taxID(sort{$a<=>$b} keys %hash_of_taxID_and_name)
							{
								my $name=$hash_of_taxID_and_name{$taxID};
								#print "Gettting full taxa for: $name\n";
								
								
								my $taxa="";
								my @arr_taxa;
								unshift(@arr_taxa,$name);
								
								my $parent_taxID=$hash_of_taxID_and_parentTaxID{$taxID};
								if($parent_taxID and $parent_taxID=~/\S+/ and $hash_of_already_created_taxaID_and_taxa{$parent_taxID})
									{
										my $e_taxa=$hash_of_already_created_taxaID_and_taxa{$parent_taxID};
										$taxa=$e_taxa.";".$name;
									}
								elsif($parent_taxID and $parent_taxID=~/\S+/)
									{	
										while($parent_taxID >1)
											{
												my $p_name=$hash_of_taxID_and_name{$parent_taxID};
												if($p_name and $p_name=~/\S+/ and $p_name!~/ group/ and $p_name!~/cellular organisms/)
													{				
														unshift(@arr_taxa,$p_name);
													}	
												
												if($hash_of_taxID_and_parentTaxID{$parent_taxID} and $hash_of_taxID_and_parentTaxID{$parent_taxID} >1)
													{
														$parent_taxID=$hash_of_taxID_and_parentTaxID{$parent_taxID};
													}
												else{
														last;
													}	
											}
											
										$taxa=join(";",@arr_taxa);	
									}
								
								#-- add a semicolon at the end
								$taxa=$taxa.";";
								
								$taxa=~s/;;/;/g;
								
								$hash_of_already_created_taxaID_and_taxa{$taxID}=$taxa;						
								print TAB "$taxID\t$taxa\n";
									
							}
						close(TAB);	
						#-------------------------------------------------------------
					}
				
						
			}
		return 1;
	}


sub identify_closest_taxa_and_get_ftps()
	{
		my($row,$hash_of_target_taxas,$hash_of_matched_taxas,$hash_of_qry_taxa_and_ftps)=@_;
		
		my @arr_r1;
		if($row=~/;|\s+/)
			{		
				@arr_r1=split(';|\s+',$row);
			}
		else{
				push(@arr_r1,$row);
			}
		
			
			
		my @arr_qry_taxa;
		foreach my $name(@arr_r1)
			{
				if($name and $name=~/\S+/ and length($name)>2)
					{
						$name=~s/\[//g;
						$name=~s/\]//g;
						
						$name=~s/\(//g;
						$name=~s/\)//g;
						
						$name=~s/\'//g;
						$name=~s/\"//g;
						
						push(@arr_qry_taxa,$name);
					}
			}
		#-- also split the last element and push into the array
		if($arr_qry_taxa[$#arr_qry_taxa]=~/\s/)
			{
				my @arr_t2=split('\s+',$arr_qry_taxa[$#arr_qry_taxa]);
				foreach my $name(@arr_t2)
					{
						if($name and $name=~/\S+/ and length($name)>2)
							{
								$name=~s/\[//g;
								$name=~s/\]//g;
								
								$name=~s/\(//g;
								$name=~s/\)//g;
								
								$name=~s/\'//g;
								$name=~s/\"//g;
								
								#push(@arr_qry_taxa,$name);
							}
					} 
			}
			
		#---- now check and score the target taxas ---
		my %hash_of_target_taxa_and_score;
		my $found=0;
		
		#--- check if multiple taxa level provided or a single one
		if($arr_qry_taxa[0] and $arr_qry_taxa[0]=~/\S+/)
			{
				#print "A:$arr_qry_taxa[0]\n";
				foreach my $taxa(sort keys %{$hash_of_target_taxas})
					{	
						
						my @arr_t1=split(';',$taxa);
						
						my %hash_of_current_taxa_element;
						foreach my $name(@arr_t1)
							{
								if($name and $name=~/\S+/ and length($name)>2)
									{
										$hash_of_current_taxa_element{$name}=1;
									}
							}
						
						#-- now score the current taxa
						#foreach my $name(@arr_qry_taxa)
						
						for(my $i=$#arr_qry_taxa; $i>=0; $i--)
							{
								my $name=$arr_qry_taxa[$i];
								
							
								if($hash_of_current_taxa_element{$name})
									{								
										if($hash_of_target_taxa_and_score{$taxa})
											{
												$hash_of_target_taxa_and_score{$taxa}=$hash_of_target_taxa_and_score{$taxa} + ($i+1);
											}
										else{
												$hash_of_target_taxa_and_score{$taxa}=$i+1;
											}
										$found++;		
									}
							}
						#print "\$found=$found\n";
					}
			}	
		else{
				return 0;			
			}	
			
		#-- if nothing matched then randomly pick any taxa
		if($found==0)
			{
				
				return 0;
			}
		
		#--- now get the ftps from all highest scoring taxa
		my $highest_score=0;
		foreach my $taxa(sort{$hash_of_target_taxa_and_score{$b}<=>$hash_of_target_taxa_and_score{$a}} keys %hash_of_target_taxa_and_score)
			{
				my $score=$hash_of_target_taxa_and_score{$taxa};
				if($score >=$highest_score)
					{
						#--
						$hash_of_matched_taxas->{$taxa}=$score;
						#-- get the list of ftps --
						my $ftps=$hash_of_target_taxas->{$taxa};
						
						my @arr_ftp;#=split(';',$ftps);
						if($ftps=~/;/)
							{
								@arr_ftp=split(';',$ftps);
							}
						else{
								push(@arr_ftp,$ftps);
							}	
						foreach my $ftp(@arr_ftp)
							{
								$hash_of_qry_taxa_and_ftps->{$row}->{$ftp}=$taxa;
							}
						#-- set the highest score
						$highest_score=$score;
					}
				if($score<$highest_score){last;}
			}
		
		return $highest_score;
	}
	
sub get_sequence_using_accessions()
	{
		my($pm,$acc_list)=@_;
		
		my @arr_accs=`cat $acc_list >&1`;
		
		#--- now download the associated protein accessions
		my @arr_commands;
		
		my $search_str_count=0;
		foreach my $search_str(@arr_accs)
			{
				chomp $search_str; $search_str=~s/\r//; $search_str=~s/\.\d$//;
				
				my $cmd1=qq~efetch -db nuccore -format fasta -id $search_str | gzip >$temp_dir/$search_str.fna.gz~;
				
				#my $a1=qq~awk '{if(\$0 ~ /\\S/){print "$search_str\\t$tab" \$0}}' ~;				
				#my $cmd2=qq~efetch -db taxonomy -id $search_str -format xml | xtract -pattern Taxon -element TaxId | awk '{if(\$0~/\S/){print $search_str"\t"\$0}else{print $search_str"\t0"}}' >$temp_dir/$search_str.fna.gz.acc_and_taxID.txt~;
				
				my $cmd2=qq~efetch -db taxonomy -id $search_str -format xml | xtract -pattern Taxon -element TaxId | awk '{if(\$0 > 0 ){print "$search_str	"\$0}else{print "$search_str	0"}}' >$temp_dir/$search_str.fna.gz.acc_and_taxID.txt~;
				
				push(@arr_commands,"$cmd1 && $cmd2");
				#push(@arr_commands2,$cmd2);
				
				$search_str_count++;
			}
								
		#---
		print "\n";	
		
		my $remaining_cmds=$search_str_count;
		foreach my $cmd(@arr_commands)
			{	
				--$remaining_cmds;
					
				print "\rRemaining: $remaining_cmds/$search_str_count \t\t\t";
				select(undef,undef,undef,0.10);
				$pm->start and next;
				
				system("$cmd");
				
				my @arr_c1=split('>',$cmd);
				my $acc_and_taxID_file=$arr_c1[$#arr_c1]; chomp $acc_and_taxID_file; $acc_and_taxID_file=~s/\r//; $acc_and_taxID_file=~s/\s+//g;
				
				#--- check the size
				my $file_size= -s $acc_and_taxID_file;
				if($file_size<1)
					{
						#print "Empty file: $acc_and_taxID_file\n";
						if($acc_and_taxID_file=~/$temp_dir\/(\S+)\.fna\.gz/)
							{
								my $acc=$1;
								system("echo '$acc	0' >$acc_and_taxID_file");
							}
						
					}		
				
				$pm->finish;
			}
		$pm->wait_all_children;
		
		print "\n";	
		
		
		#------- some accessions returns no TaxID and xtract terminates without returning anything, check and fix them
		 
		return 1;
	}	

	
	
sub get_unique_id()
	{
		my $letters="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
		my @arr_letters=split('',$letters);		
		my $unique_id = time();	   
		for(my $i = 0; $i < int(rand(50)); $i++)
			 {
				$unique_id .= $arr_letters[int(rand($#arr_letters))];
			 }
		return $unique_id;	 
	}

	
