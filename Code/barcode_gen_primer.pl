#########################################
# This file is part of barCoder.
#
# barCoder is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# barCoder is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with barCoder.  If not, see <https://www.gnu.org/licenses/>.
#
#########################################
# Subroutine genPrimer: produces a primer sequence meeting the parameters given as 
#  input. There are two categories of these parameters, as well as some general
#  inputs. These inputs are all required:
# Inputs:
#   (1) General parameters
#      - current target genome name: needed for labeling outputs
#      - current primer type: used primarily for keeping track of which 
#		parameters are which in the final output; also in some 
#		differences in PCR requirements for primers/probes
#      - v: toggle for verbose console output
#   (2) Blast parameters
#      - blast threshold: the hit score cutoff for when a primer matches 
#		somewhere on the genome with too high of a score
#   (3) PCR parameters
#      - primer length (min/max): PCR parameters
#      - primer target Tm (min/max): PCR parameters
#      - primer GC content (min/max): PCR parameters
#      - # of consecutive A/T/C's allowed: probe PCR parameter
#      - # of consecutive G's allowed: probe PCR parameter
# Outputs:
#   - sequence of randomly generated primer that passes all tests
#########################################
sub genPrimer {
#########################################
# Interpret input vars
#########################################
#General Parameters
my $curRep=@_[0];		#counter for generating multiple alternative
				# barcodes
my $curGenome=@_[1];		#name of the file with info on the current genome
my $curPrimerType=@_[2];	#project primer, target primer, or target probe
my $v=@_[3]; 			#verbose output off or on
my $log=@_[4];			#write to log file on or off
#Blast Parameters
my $blastDB=@_[5];		#location of NCBI db
my $blastThresh=@_[6];		#max allowed blast hit
my $blastCheckDb=@_[7];         #whether to check NCBI db
#Primer Parameters
my $lengthMin=@_[8];		#minimum allowed primer length 
my $lengthMax=@_[9];		#maximum allowed primer length
my $TmMin=@_[10];		#minimum allowed Tm
my $TmMax=@_[11];		#maximum allowed Tm
my $GCmin=@_[12];		#minimum allowed GC-content
my $GCmax=@_[13];		#maximum allowed GC-content
my $maxATCreps=@_[14]+1;	#maximum number of A/T/C bp repeats
my $maxGreps=@_[15]+1;		#maximum number of allowed G repeats
my $minHbonds=@_[16];		#minimum number of hyrdogen bonds in stem of 
				# stem-loop to be a problem
my $minPalLength=@_[17];	#minimum palindrome length to be a problem
my $maxPalLength=@_[18];	#maximum palindrome length to be a problem
my $gapLimit=@_[19];		#largest loop length in stem-loop to be a problem
my $numMismatch=@_[20];		#maximum number of mismatches in stem of
				# stem-loop to be a problem
chomp(my $cpu_count = `grep -c -P '^processor\\s+:' /proc/cpuinfo`); # cpu count
#print $cpu_count;
#########################################
# Basic declarations, etc.
#########################################
my $projPath=getcwd();

#########################################
# Compute the viable options for # of A/T's and # of G/C's to meet length, Tm,
# and %GC-content requirements
#########################################
my @probeParams=();
for my $length ($lengthMin ... $lengthMax) {
	for my $numGC (0 ... $length) {
		$percGC=100*$numGC/$length;
		$Tm=69.4+41*($numGC-16.4)/$length;
		if (($percGC>=$GCmin) && ($percGC<=$GCmax) &&	
		   ($Tm>=$TmMin) && ($Tm<=$TmMax)) {
			push(@probeParams,$length-$numGC,$numGC);
		}
	}
}

#########################################
# Main loop: cycles through the various requirements, fixing where possible,
# or starting with a new random sequence if necessary.
#########################################
#update log and user on progress
update("\nSearching for $curPrimerType for target genome file $curGenome:\n",
	$v,$log);
#while loop makes sure the code runs until a valid sequence is generated;
#if statement allows it to break the loop and try again
my $passAll=0;
while ($passAll==0) { 
	#########################################
	# Randomly generate a primer by (1) randomly selecting the 
	# AT vs GC content and length based on precomputed viable sets 
	# stored in @probe_params, and (2) generating a random sequence 
	# w/ the right number of ATs/GCs
	#########################################
	# Randomly select AT vs GC content from precomputed table
	my $nOptions = @probeParams/2;
	my $randi=int(rand($nOptions));
	# Call subroutine to generate a random sequence w/ the right number 
	# of ATs/GCs
	my $nAT=$probeParams[2*$randi];
	my $nGC=$probeParams[2*$randi+1];
	my $randSeq=genRand($nAT,$nGC);
	#update log and user on progress
	update("\tNew sequence generated\n"."Current sequence = ".$randSeq.
		"\n",$v,$log);
	#########################################
	# Check if sequence matches requirements, specifically:
	#   - More C's than G's
	#   - No G on 5' end
	#   - <= 3 consecutive Gs, <=4 A/T/Cs
	#   - Does not form primer dimers with other primers in the project
	#   - Does not have homology with the other primers in the project
	#   - Does not have homology with any of the project-relevant genomes
	#########################################
	#----------------------------------------
	# Make sure there are more Cs than Gs; take complement otherwise (only
	# relevant for probes, not primers)
	#----------------------------------------
	if ($curPrimerType=~m/[Pp]robe/) {
		#update log and user on progress
		update("\tChecking number of G's vs C's (probe only)...",$v,$log);
		# First, count number of each
		my $nCs=0; $nCs++ while($randSeq =~ m/C/g);
		my $nGs=0; $nGs++ while($randSeq =~ m/G/g);
		# Now, if more Gs than Cs, simply take the complement
		if ($nCs<$nGs) {
			$randSeq=seqComp($randSeq);
			#update log and user on progress
			update("\n\t\tFAIL: More Gs than Cs found, ".
				"complement taken\n".
				"Current sequence = ".$randSeq."\n",$v,$log);
		} else {
			#update log and user on progress
			update("OK\n",$v,$log);
		}
	}
	#----------------------------------------
	# While loop to check for fixable issues
	#----------------------------------------
	my $fixableIssues=1;
	while ($fixableIssues==1) {
		#----------------------------------------
		# Check for G on 5' end; shift sequence if needed (only 
		# relevant for probes, not primers)
		#----------------------------------------
		if ($curPrimerType=~m/[Pp]robe/) {
			#update log and user on progress
			update("\tChecking for G on 5' end (probe only)...",
				$v,$log);
			#split into array
			my @bps=split('',$randSeq);
			#while there is still a G at the 1st postion, 
			#shift it to the end
			if ($bps[0] eq "G") {
				shift(@bps);
				push(@bps,'G');
				#rejoin back into a string
				$randSeq=join('',@bps);
				#update log and user on progress
				update("\n\t\tFAIL: G found @ 5' end, sequence ".
					"shifted\n"."Current sequence = ".
					$randSeq."\n",$v,$log);
				#try again
				redo;
			}  else {
				#update log and user on progress
				update("OK\n",$v,$log);
			}
		} 
		#----------------------------------------
		# Check for # G/C's !=2 in last 5 bases at 3' end (only relevant  
		# for primers, not probes
		#----------------------------------------		
		else {
			#update log and user on progress
			update("\tChecking for # C/G's !=2 in last 5 bps @ 3' end ".
				"(primers only)...",$v,$log);
			#split into array
			my @bps=split('',$randSeq);
			#count G's and C's at 3' end
			my $GCcnt=0;
			for my $bp (1...5) {
				if (@bps[-$bp]=~/[GgCc]/) { $GCcnt++; }
			}
			#while there is still a G at the 1st postion, 
			#shift it to the end
			if ($GCcnt != 2) { 
				#shuffle sequence
				$randSeq=shuffleSeq($randSeq);
				#update log and user on progress
				update("\n\t\tFAIL: ".$GCcnt." C/G's found, ".
					"sequence shuffled\n".
					"Current sequence = ".
					$randSeq."\n",$v,$log);
				#try again
				redo;
			} else {
				#update log and user on progress
				update("OK\n",$v,$log);
			}
		}
		#----------------------------------------
		# Check for repeats (typically <= 3 Gs, <= 4 A/T/C's)
		# and shuffle if needed
		#----------------------------------------
		#update log and user on progress
		update("\tChecking for bp repeats...",$v,$log);
		#if there are still repeats, flag which bp it was
		my $bpRepeatFlag=0;
		if ($randSeq=~m/A{$maxATCreps,}/) { $bpRepeatFlag='A'; }
		elsif ($randSeq=~m/T{$maxATCreps,}/) { $bpRepeatFlag='T'; }
		elsif ($randSeq=~m/C{$maxATCreps,}/) { $bpRepeatFlag='C'; }
		elsif ($randSeq=~m/G{$maxGreps,}/) {$bpRepeatFlag='G'; } 
		#if the flag is thrown shuffle and report which bp it was
                if ($bpRepeatFlag=~m/\D/) {
			#shuffle sequence
			$randSeq=shuffleSeq($randSeq);
			#update log and user on progress
			update("\n\t\tFAIL: Too many consecutive ".$bpRepeatFlag.
				"'s found, sequence shuffled\n"."Current ".
				"sequence = ".$randSeq."\n",$v,$log);
			#try again
			redo;
		} else {
			#update log and user on progress
			update("OK\n",$v,$log);
		}
		#----------------------------------------
		# Check for start codons
		#----------------------------------------		
		#update log and user on progress
		update("\tChecking for start codons...",$v,$log);
		#check for start codons (including reverse complements)
#		if ($randSeq=~/(ATG|CTG|TTG|CAT|CAG|CAA)/i) {
		if ($randSeq=~/(ATG|CAT)/i) {
			#shuffle sequence
			$randSeq=shuffleSeq($randSeq);
			#update log and user on progress
			update("\n\t\tFAIL: Start codon $1 found, sequence ". 
				"shuffled\n"."Current sequence = ".$randSeq.
				"\n",$v,$log);
			#try again
			redo;
		} else {
			#update log and user on progress
			update("OK\n",$v,$log);
		}
		#----------------------------------------
		# Check for obvious stem-loop structures
		#----------------------------------------
		#update log and user on progress
		update("\tChecking for stem-loops...",$v,$log);
		#If stem-loop detected, reshuffle
		if (stemLoopCheck($randSeq,$minHbonds,$minPalLength,$maxPalLength,
				$gapLimit,$numMismatch)) {
			#shuffle sequence
			$randSeq=shuffleSeq($randSeq);
			#update log and user on progress, as necessary
			update("\n\t\tFAIL: Stem-loop structure detected, ".
				"sequence shuffled\n"."Current sequence = ".
				$randSeq."\n",$v,$log);
			#try again
			redo;
		} else {
			#update log and user on progress
			update("OK\n",$v,$log);
		}
		#----------------------------------------
		# Done with PCR checks 
		#----------------------------------------		
		$fixableIssues=0;
	}
	#----------------------------------------
	# Check for blast hits against:
	#   - list of previous primers 
	#   - list of primers generated so far for the current project
	#   - target genomes
	#   - other genomes
	#----------------------------------------
	#Put random sequence in seq object form
	my $randSeqObj=Bio::Seq->new(-id=>'rand_primer',-seq => $randSeq);
	#Check current primer against primers from previous projects
	my $subPath = '/allPrevPrimers.fa';
	my $pass=0;
	$maxScore=blastFile($subPath,$blastThresh,$randSeqObj, $cpu_count);
	#update log and user on progress
	update("\tBlasting vs previous primers file: ".
		"allPrevPrimers.fa...",$v,$log);
	#if sequence passes being blasted against list of previous projects'
	# primers, check list of current project's primers
	if ($maxScore<=$blastThresh) {
		#update log and user on progress
		update("OK\n"."\tBlasting vs ".
			"current primers file: primerList.fa...",$v,$log);
		#Check current primer against the current list of primers for the
		# current project
		$subPath = '/primerList.fa';
		$maxScore=blastFile($subPath,$blastThresh,$randSeqObj, $cpu_count);
	} 
	#if sequence passes being blasted against list of current project's
	# primers, check against target genomes
	if ($maxScore<=$blastThresh) { 
		#update log and user on progress
		update("OK\n",$v,$log);
		#Next check against the target genomes for the project
		$maxScore=blastGenomeDir(
			'targetGenomes',$blastThresh,$randSeqObj,$v,$log, $cpu_count);
	}
	#if sequence passes being blasted against target genomes, check against
	# other genomes
	if ($maxScore<=$blastThresh) {
		#blast primer against other genomes
		$maxScore=blastGenomeDir(
			'otherGenomes',$blastThresh,$randSeqObj,$v,$log, $cpu_count);
	}
	#if sequence passes being blasted against other genomes, check against
	# NCBI db
	if ($maxScore<=$blastThresh && $blastCheckDb) {
		#update log and user on progress
		update("\tBlasting vs NCBI database...",$v,$log);
		#set up blast factory for NCBI db
		my @blastDBparams=(-db_name=>$blastDB);
		my $blastDBfactory=
			Bio::Tools::Run::StandAloneBlastPlus->new(@blastDBparams);
		#run blastall against the db
		
		if (length $randSeqObj < 20) { update("\nWarning: BLAST shorter that 20 bases\n",$v,$log);}
		my $blastDBreport=$blastDBfactory->blastn(-query=>$randSeqObj, 
							-method_args => [ -task => 'blastn-short', -num_threads => "\'".$cpu_count."\'"]);
		#grab the first result
		#my $topBlastDBresult=$blastDBreport->next_result;
		my $topBlastDBresult=$blastDBreport->next_hit;
		#if there are no results, score = 0
		#if ($topBlastDBresult->num_hits==0) {
		if (!defined($topBlastDBresult)) {
			$maxScore=0;
		#else, grab the score of the top hit
		} else {
			#grab top hit
			my $topBlastDBhit=$topBlastDBresult;
			#compute score of the top hit
			$maxScore=$topBlastDBhit->raw_score/$randSeqObj->length;
		}
		$blastDBfactory->cleanup;
	}	
	#if sequence doesn't pass any of the blast checks, try again
	if ($maxScore>$blastThresh) {
		#update log and user on progress
		update(sprintf("\n\t\tFAIL: Blast hit > %.0f%% similarity found,".
			" trying again\n",100*$blastThresh),$v,$log);
		redo;
	} 		
	#----------------------------------------
	# Primer sequence passed all the requirements! Now need to save it to
	# the primerList.fa file.
	#----------------------------------------	
	#update log and user on progress, as necessary
	update("Sequence passed! ".$randSeq."\n",$v,$log);
	#Update primerList.fa file
	open(curPrimerList,">>primerList.fa") or die ("Failed to open: $!\n");
	print curPrimerList "> $curGenome, $curPrimerType, Rep=$curRep\n".
		$randSeq."\n";
	close(curPrimerList);
	#end loop
	$passAll=1;
	return $randSeq;
} #end while and if statements
} #end genPrimer subroutine

#########################################
# Subroutine genRand: generates a random primer sequence
# Inputs:
#   - number of As/Ts
#   - number of Gs/Cs
# Outputs:
#   - random primer sequence
#########################################
sub genRand {
	# Rename input vars
	my $nAT=$_[0];my $nGC=$_[1];
	# Declare output var
	my $randSeq=();
	# Declare temp vars
	my $i;my $rAT;my $rGC;
	# Generate a random string of ATs of length nAT
	for ($i=0;$i<$nAT;$i++) {
		$rAT=int(rand(2));
		if ($rAT==0) { $randSeq=$randSeq.'A'; }
		if ($rAT==1) { $randSeq=$randSeq.'T'; }
	}
	#generate a random string of GCs of length nGC
	for ($i=0;$i<$nGC;$i++) {
		$rGC=int(rand(2));
		if ($rGC==0) { $randSeq=$randSeq."G"; }
		if ($rGC==1) { $randSeq=$randSeq."C"; }
	}
	#shuffle the sequence 
	my $randSeqShuffled=shuffleSeq($randSeq);
	return $randSeqShuffled;
} #end genRand subroutine

#########################################
# Subroutine shuffleSeq: shuffles the order of a DNA sequence
# Inputs:
#   - sequence to be shuffled
# Outputs:
#   - shuffled sequence
#########################################
sub shuffleSeq {
	#rename input var
	my $randSeq=join('',@_);
	#split into array
	my @bps=split('',$randSeq); 
	#shuffle array
	my @bpsShuffled=shuffle(@bps); 
	#rejoin back into a string
	my $randSeqShuffle=join('',@bpsShuffled); #rejoin string
	return $randSeqShuffle;
} #end shuffleSeq subroutine

#########################################
# Subroutine seqComp: takes the complement of a sequence
# Inputs:
#   - sequence to be complemented
# Outputs:
#   - complemented sequence
#########################################
sub seqComp {
	#rename input var
	my $randSeq=join('',@_);
	#split into array
	my @seq1=split('',$randSeq);
	#declare output var 
	my @seq2=();
	#run through each base and take the complement
	foreach my $bp1 (@seq1) {
		my $bp2=0;
		if ($bp1=~m/C/g) {
			$bp2='G';
		} elsif ($bp1=~m/G/g) {
			$bp2='C';
		} elsif ($bp1=~m/A/g) {
			$bp2='T';
		} elsif ($bp1=~m/T/g) {
			$bp2='A';
		}
		#store the complement bp
		push(@seq2,$bp2);
	}
	#rejoin back into a string
	$randSeq=join('',@seq2);
	return $randSeq;
} #end seqComp subroutine

#########################################
# Subroutine stemLoopCheck: looks for potential stem-loop structures in a primer
# Inputs:
#   - sequence of current primer
# Outputs:
#   - 1 if stem-loop structure detected
#   - 0 otherwise
#########################################
sub stemLoopCheck {
	#interpret inputs
	my $randSeq=@_[0];
	my $minHbonds=@_[1];
	my $minPalLength=@_[2];
	my $maxPalLength=@_[3];
	my $gapLimit=@_[4];
	my $numMismatch=@_[5];
	#initialize variables
	my $stemLoopDetected=0;
	my $randSeqObj=Bio::Seq->new(-id=>'rand_primer',-seq => $randSeq);
	#initialize EMBOSS palindrome factory
	my $EMBOSSfactory=Bio::Factory::EMBOSS->new();
	my $palindrome=$EMBOSSfactory->program('palindrome');
	#run EMBOSS palindrome algorithm
	my $palResults=$palindrome->run({-sequence=>$randSeqObj,
					 -minpallen=>$minPalLength,
					 -gaplimit=>$gapLimit,
					 -maxpallen=>floor($randSeqObj->length/2),
					 -nummismatches=>$numMismatch,
					 -outfile=>'palResults.out'});
	#interpret results
	my $parsedResults=new Bio::Tools::EMBOSS::Palindrome(
		-file => 'palResults.out');
	while (my $seq=$parsedResults->next_seq) {
		for my $feat ($seq->get_SeqFeatures) {
			#grab key palindrome features
			my $palTopStart=$feat->start;
			my $palTopEnd=$feat->end;
			my $palBotStart=$feat->hstart;
			my $palBotEnd=$feat->hend;
			my $palSeq=$randSeqObj->subseq($palTopStart,$palTopEnd);
			my $palGap=$palBotStart-$palTopEnd;
			#compute number of H bonds
			my $palGCs=0; while($palSeq=~m/[GC]/g) { $palGCs++ };
			my $palATs=$feat->length-$palGCs;
			my $palHbonds=3*$palGCs+2*$palATs;
			#check for stem-loop
			if ($palHbonds>$minHbonds && $palGap<$gapLimit) { 
				$stemLoopDetected=1; 
				last;
			}
		}
	}
	unlink('palResults.out');
	return $stemLoopDetected;
} #end primerDimer subroutine

#########################################
# Subroutine blastGenome: checks a primer for blast hits against any genomes 
#  contained in a specific folder.
# Inputs:
#   - subdir: name of folder containing the files of genomes to be checked
#   - threshold: value for the maximum hit score that's acceptable. This allows
# 	the script to abort upon failure rather than continuing to blast every #	genome every time.
#   - randSeqObj: sequence to blast
#   - v: toggle for verbose output
# Outputs:
#   - returns 1 if sequence passes, 0 if sequence fails
#########################################
sub blastGenomeDir {
	#interpret inputs
	my $subdir=@_[0];
	my $threshold=@_[1];
	my $randSeqObj=@_[2];
	my $v=@_[3];
	my $log=@_[4];
	my $cpu_count=@_[5];
	#declare key variables
	my $curScore=0;
	#change directory grab filenames
	my $projPath=getcwd();
	chdir("$projPath/$subdir") or die ("Failed to change directories: $!\n");
	my @genomeFiles=<*>;
	chdir($projPath);
	#go through directory and check for hits
	foreach my $genomeFile (@genomeFiles) {
		#don't blast deleted files
		if ($genomeFile=~m/~$/) { next; }
		#compute path to current file
		my $subPath="/$subdir/$genomeFile";
		#update log and user on progress
		update("\tBlasting vs $subdir genome file: $genomeFile...",
			$v,$log);
		#cycle through all the sequences in a file
		$maxScore=blastFile($subPath,$threshold,$randSeqObj, $cpu_count);
		#if it's a fail, break all the way out
		if ($maxScore>$threshold) { 
			last; 	#abort current folder
		} else {
			#update log and user on progress
			update("OK\n",$v,$log);
		}
	}
	#go back to original directory
	chdir($projPath);
	#update output
	return $maxScore;
} #end blastGenomeDir subroutine

#########################################
# Subroutine blastFile: Blasts a primer against all sequences in a file.
# Inputs:
#   - path: path to file containing the sequences to be scanned
#   - threshold: value for the maximum hit score that's acceptable. This allows
# 	the script to abort upon failure rather than continuing to blast every #	genome every time.
#   - randSeqObj: sequence to blast
#   - v: toggle for verbose output
# Outputs:
#   - returns 1 if sequence passes, 0 if sequence fails
#########################################
sub blastFile {
	#interpret inputs
	my $subPath=@_[0];
	my $threshold=@_[1];
	my $randSeqObj=@_[2];
	my $cpu_count = @_[3];
	#declare key variables
	my $projPath=getcwd();
	#load sequence file
	my $seqFile=Bio::SeqIO->new(-id=>'curSeqFile',-file=>$projPath.$subPath);
	# Set blast parameters
	#my @blastParams=(program=>'blastn');
	my $blastFactory=Bio::Tools::Run::StandAloneBlastPlus->new();
	#go through each sequence
	my $moreSeqs=1;
	my $curScore=0;
	my $maxScore=0;
	while ($moreSeqs==1) {
		#grab next sequence
		my $curSeq=$seqFile->next_seq;
		#make sure there's still a sequence, otherwise end loop
		if ($curSeq=~m/^$/) { last; }
		#blast primer sequence against current file sequence
		if (length $randSeqObj < 20) { update("\nWarning: BLAST shorter that 20 bases\n",$v,$log);}

		my $curReport=$blastFactory->bl2seq(-method=>'blastn', -query=>$randSeqObj, -subject=>$curSeq,
						    -method_args => [ -task => 'blastn-short', -num_threads => "\'".$cpu_count."\'"]);
		#process blast report
		#my $curResult=$curReport->next_iteration;
		#make sure there are hits; if not, go to next sequence in file
		if ($curReport->num_hits==0) { redo; }
		#load first hit and compute the score
		my $curHit=$curReport->next_hit;
		$curScore=$curHit->raw_score/$randSeqObj->length;
		#if hit is over threshold, abort
		if ($curScore > $threshold) {
			$maxScore=$curScore;
			last;		#abort current file
		#else record if current score is the biggest yet (for reporting
		# purposes only)
		} elsif ($curScore>$maxScore) {
			$maxScore=$curScore;
		} 
	}
	$blastFactory->cleanup;  #Cleanup	temporary files
	return $maxScore;
} #end blastFile subroutine

#########################################
# Subroutine update: Updates log file and user on progress, depending on 
#  options chosen by the user via the parameters.txt file.
# Inputs:
#   - update: string containing the update text
#   - v: whether or not to output updates to screen
#   - log: whether or not to output updates to log file
# Outputs:
#   - none
#########################################
sub update {	
	#interpret inputs
	my $update=@_[0];	#string containing update text
	my $v=@_[1];		#output to screen flag
	my $log=@_[2];		#output to log file flag
	#update log and user on progress, as necessary
	if ($log==1) { print logFile $update; }
	if ($v==1) { print $update; }
} #end update subroutine

1;	
