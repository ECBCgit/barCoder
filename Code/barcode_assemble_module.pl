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
# Subroutine addSpacer: generates spacer sequence to fill in barcode module. The
#  spacer sequence is chosen to match the GC content of the organism and to create
#  no obvious stem-loop structures.
# Inputs:
#   - genome sequence filename
#   - sequence of the 2 primers and probe for the current barcode. These sequences
#	are in a string, not a seq object.
#   - general parameters
#   - stem-loop detection parameters
# Outputs:
#   - write barcode sequence to file barcodeList.fa
#########################################
sub addSpacer {
	#interpret inputs
	my $curRep=@_[0];	#counter for generating multiple alternative
				# barcodes
	my $genome=@_[1];	#file containing current target genome
	my $primer1=@_[2];	#sequence of primer1 for current barcode
	my $primer2=@_[3];	#sequence of primer2 for current barcode
	my $probe=@_[4];	#sequence of probe for current barcode
	my $v=@_[5];		#flag to output progress to screen
	my $log=@_[6];		#flag to output progress to log file
	my $minHbonds=@_[7];	#minimum number of Hyrdogen bonds in stem of 
				# stem-loop to be a problem
	my $minPalLength=@_[8];	#minimum palindrome length to be a problem
	my $maxPalLength=@_[9];	#maximum palindrome length to be a problem
	my $gapLimit=@_[10];	#largest loop length in stem-loop to be a problem
	my $numMismatch=@_[11];	#maximum number of mismatches in stem of
				# stem-loop to be a problem
	my $primerGap=@_[12];	#number of base pairs separating primers (includes
				# the length of the probe)
	#------------------------------
	# Match spacer GC content to genome GC content
	#------------------------------
	#update log and user on progress
	update("\nAll primers generated, adding spacers: \n".
		"\tComputing spacer GC content to match genome...",$v,$log);
	#compute the GC content of the target genome and 3 primers
	my $percentGenomeGC=getFileGCcontent($genome);
	my $GCprimer1=0; while($primer1=~m/[GC]/g) { $GCprimer1++ };
	my $GCprimer2=0; while($primer2=~m/[GC]/g) { $GCprimer2++ };
	my $GCprobe=0; while($probe=~m/[GC]/g) { $GCprobe++ };
	#compute number of GCs and total length needed for the spacers
	my $totalBClength=$primerGap+length($primer1)+length($primer2);
	my $totalGCsNeeded=floor($percentGenomeGC*$totalBClength);
	my $spacerGCsNeeded=$totalGCsNeeded-$GCprimer1-$GCprimer2-$GCprobe;
	my $totSpacerLength=$primerGap-length($probe);
	my $GCsPerProbe=floor($spacerGCsNeeded/2);
	my $ATsPerProbe=floor(($totSpacerLength-$spacerGCsNeeded)/2);
	#update log and user on progress
	update("done\n",$v,$log);
	#------------------------------
	# Generate random spacer sequences and check for stem-loops and 
	# start codons
	#------------------------------
	my $barcodePass=0;
	my $spacer1;
	my $spacer2;
	my $barcode;
	while ($barcodePass==0) {
		#set the flag to pass unless it fails below
		$barcodePass=1;
		#generate spacer sequence
		$spacer1=genRand($ATsPerProbe,$GCsPerProbe);
		$spacer2=genRand($ATsPerProbe,$GCsPerProbe);
		#take reverse complement of primer 2
		my $primer2rc=reverse($primer2);
		$primer2rc =~ tr/ATCGatcg/TAGCtagc/;		
		#assemble barcode into a single sequence
		$barcode=$primer1.$spacer1.$probe.$spacer2.$primer2rc;
		#update log and user on progress
		update("Current barcode = $barcode\n",$v,$log);
		#check for stem-loops
		$stemLoopDetected=stemLoopCheck($barcode,$minHbonds,
			$minPalLength,$maxPalLength,$gapLimit,$numMismatch);
		#update log and user on progress
		update("\tChecking for stem-loops...",$v,$log);
		if ($stemLoopDetected) {
			#update log and user on progress
			update("\n\t\tFAIL: stem-loop found, trying new spacers\n",
				$v,$log);
			#set the pass flag to zero
			$barcodePass=0;
		} 
		#if no stem-loop detected
		if ($barcodePass==1) {
			#update log and user on progress
			update("OK\n"."\tChecking for start codons...",$v,$log);
			#check for start codons
#			if ($barcode=~/(ATG|CTG|TTG|CAT|CAG|CAA)/i) {
			if ($barcode=~/(ATG|CAT)/i) {
				#update log and user on progress
				update("\n\t\tFAIL: start codon $1 found, ".
					"trying new spacers\n",$v,$log);
				#set the pass flag to zero
				$barcodePass=0;
			}
		}
	}
	#------------------------------
	# Barcode complete; wrap things up
	#------------------------------	
	#update log and user on progress
	update("OK\n",$v,$log);
	#now write entire barcode to file and return it
	#Update barcodeList.fa file
	open(curBarcodeList,">>barcodeList.fa") or die ("Failed to open: $!\n");
	print curBarcodeList "> $genome, $curPrimerType, Rep=$curRep\n".
		$barcode."\n";
	close(curBarcodeList);
	#update log and user on progress,
	update("\nComplete barcode generated for $genome! \n\n".
		"Complete barcode sequence = \n\n".$barcode."\n",$v,$log);
	#return barcode sequence
	return $barcode;
} #end addSpacer subroutine

#########################################
# Subroutine getGCcontent: grabs the GC content from a genome sequence file
# Inputs:
#   - genome sequence file
# Outputs:
#   - proportion of sequence composed of G/C's (0-1)
#########################################
sub getFileGCcontent {
	#interpret inputs
	my $genomeFile=@_[0];
	#open sequence file
	my $genomeSeq=Bio::SeqIO->new(-file=>'targetGenomes/'.$genomeFile);
	#go through each entry in the file and count the total GC and total length
	my $moreSeqs=1;
	my $totGCs=0;
	my $totLength=0;
	while ($moreSeqs==1) {
		#grab next sequence
		my $curSeqObj=$genomeSeq->next_seq;
		#make sure there's still a sequence, otherwise end loop
		if ($curSeqObj=~m/^$/) { last; }
		#count the G/C's
		my $curSeqSeq=$curSeqObj->seq;
		my $curLength=$curSeqObj->length;
		my $curGCs=0; while($curSeqSeq=~m/[GC]/g) { $curGCs++ };
		$totGCs=$totGCs+$curGCs;
		$totLength=$totLength+$curLength;
	}
	#compute the GC content from the total GC count and the total length
	return $totGCs/$totLength;
} #end getGCcontent subroutine

1;
