#!/usr/bin/perl
#########################################
# This file is part of Barcoder.
#
# Barcoder is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# Barcoder is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Barcoder.  If not, see <https://www.gnu.org/licenses/>.
#
#########################################
use strict;
use warnings;
use Cwd;
use POSIX;
use List::Util qw(min max shuffle);
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Bio::Index::Fasta;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::Tools::EMBOSS::Palindrome;
use Bio::Tools::Run::EMBOSSApplication;
use Bio::Tools::GFF;
use Bio::Factory::EMBOSS;
use Bio::AlignIO;
#########################################
# Pull in other perl scripts
#########################################
require 'Code/barcode_gen_primer.pl';
require 'Code/barcode_extract_params.pl';
require 'Code/barcode_assemble_module.pl';
#########################################
# Load input information
#########################################
# Change directory to correct project folder
my $projectFolder=$ARGV[0];
my $nBarcodes=$ARGV[1];
my $targetFlag=$ARGV[2];
chdir("Projects/".$projectFolder);
my $projPath=getcwd();
#----------------------------------------
# Load parameters from parameter file
#----------------------------------------
my @params=extractParams();
# Declare parameter groups
my @initParams=();
my @genParams=();
my @blastParams=();
my @primer1Params=();
my @primer2Params=();
my @probeParams=();
my @spacerParams=();
my @targetParams=();
# Fill in each parameter group
my $groupCnt=0; my $paramCnt=0;
foreach my $param (@params) {
	if ($param=~/x/) { $groupCnt++; $paramCnt=0 }
	elsif ($groupCnt==0) { $initParams[$paramCnt]=$param; $paramCnt++; }
	elsif ($groupCnt==1) { $genParams[$paramCnt]=$param; $paramCnt++; }
	elsif ($groupCnt==2) { $blastParams[$paramCnt]=$param; $paramCnt++; }
	elsif ($groupCnt==3) { $primer1Params[$paramCnt]=$param; $paramCnt++; }
	elsif ($groupCnt==4) { $primer2Params[$paramCnt]=$param; $paramCnt++; }
	elsif ($groupCnt==5) { $probeParams[$paramCnt]=$param; $paramCnt++; }
	elsif ($groupCnt==6) { $spacerParams[$paramCnt]=$param; $paramCnt++; }
	elsif ($groupCnt==7) { $targetParams[$paramCnt]=$param; $paramCnt++; }
}
#----------------------------------------
# Set up logs
#----------------------------------------
chdir($projPath."/Logs") or die "Failed to open directory Logs: $!";
# If instructed to, remove old log files
my $clearLogsFlag=$initParams[0];
if ($clearLogsFlag==1) { 
	foreach my $logFile (<*>) 
		{ unlink($logFile); 
	} 
}
# Grab system time
(my $sec,my $min,my $hour,my $mday,my $mon,my $year,my $wday,my $yday,
	my $isdst)=localtime(time);
# Format timestamp
my $timeStamp=sprintf "%4d-%02d-%02d %02d:%02d:%02d",
		$year+1900,$mon+1,$mday,$hour,$min,$sec;
# Open log file
open(logFile,">".$timeStamp.".log") or die "Failed to open results file: $!";
chdir($projPath);
#----------------------------------------
# Load any primers from previous projects to exclude
#----------------------------------------
# First clear out any old file contents of the master list
open(allPrevPrimers,">allPrevPrimers.fa") 
	or die "Failed to open results file: $!";
print allPrevPrimers "";
close(allPrevPrimers);
# Now reopen so that the new information can be added
open(allPrevPrimers,">>allPrevPrimers.fa") 
	or die "Failed to open results file: $!";
# Switch to directory containing primers to avoid
chdir($projPath."/previousPrimers") or die "Failed to open directory Logs: $!";
# Run through all files containing lists of previous primers
foreach my $file (<*>) {
	open(prevPrimerList,$file) or die "Failed to open file: $!";
	while (<prevPrimerList>) {
		print allPrevPrimers $_;
	}
	close(prevPrimerList);
}
close(allPrevPrimers);
chdir($projPath);
#########################################
# Make sure output files exist and are cleared (as appropriate)
#########################################
#Make sure primerList.fa exists, and clear if it does
open(primerList, ">primerList.fa") or die ("Failed to load file: $!\n");
print primerList "";
close(primerList);
#Make sure primerList.fa exists, and clear if it does
open(curBarcodeList,">barcodeList.fa") or die ("Failed to open: $!\n");
print curBarcodeList "";
close(curBarcodeList);
#########################################
# Generate barcode modules. Each module includes 3 primers (project 
# primer, target primer, and target probe). The 3 primers are connected by spacer
# sequences to make up the entire barcode module. The primers meet a set of
# requirements to maximize the chance that (1) the primers will work well and (2) 
# the primers are appropriately unique. The resulting primers are stored in 
# primerList.fa and the full barcodes are stored in barcodeList.fa
#########################################
my $dir=getcwd();
chdir('targetGenomes');
my @genomes=<*>;
chdir($dir);
# Repeat so that you get a number of barcodes = nBarcodes for each target
# genome
for (my $j=1; $j<=$nBarcodes; $j++) {
	#Generate the project primer
#	my $genomeFile='all genomes in project';
#	my $primerType='Project Primer';
#	my $projPrimerSeq=genPrimer($j,$genomeFile,$primerType,@genParams,
#		@blastParams,@projPrimerParams);
	#go through each target genome
	foreach my $genome (@genomes) {
		#skip files ending in ~
		if ($genome=~/\w$/) {
			# Generate target primer 1	
			my $primerType='Target Primer 1';
			my $primer1Seq=genPrimer($j,$genome,$primerType,
				@genParams,@blastParams,@primer1Params);	
			# Generate target primer 2	
			$primerType='Target Primer 2';
			my $primer2Seq=genPrimer($j,$genome,$primerType,
				@genParams,@blastParams,@primer2Params);	
			# Generate target probe
			$primerType='Target Probe';
			my $probeSeq=genPrimer($j,$genome,$primerType,
				@genParams,@blastParams,@probeParams);
			# Generate spacer sequence and write barcode to file
			addSpacer($j,$genome,$primer1Seq,
				$primer2Seq,$probeSeq,@genParams,
				@spacerParams);
		}
	}
}

#########################################
# Wrap everything up
#########################################
print logFile "\nScript end.\n\n";
print "\nScript end.\n\n";
close(logFile);


