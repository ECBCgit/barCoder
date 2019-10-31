#!/usr/local/bin/perl
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
use strict;
use warnings;

print "";

sub extractParams {

my @initParams=();
my @genParams=();
my @blastParams=();
my @primer1Params=();
my @primer2Params=();
my @probeParams=();
my @spacerParams=();
my @targetParams=();

open(paramFile,"parameters.txt") or die "Failed to open parameters file: $!";
while (<paramFile>) {
	my $curLine=$_;
	#First grab initialization parameters
	if ($curLine =~ m/initialize_(\S*) = (\d*)/) {
		my $hit=$2;
		if ($1=~/clearLogs/) { $initParams[0]=$hit; }
		elsif ($1=~/EMBOSSpath/) { system("export PATH=\$PATH:$hit"); }
	#Now grab general parameters
	} elsif ($curLine =~ m/general_(\w*) = (\d*)/) {
		my $hit=$2;
		if ($1=~/verbose/) { $genParams[0]=$hit; }
		elsif ($1=~/logs/) { $genParams[1]=$hit; }	
	#Now grab blast parameters
	} elsif ($curLine =~ m/blast_(\w*) = ([\S]*)/) {
		my $hit=$2;
#		if ($1=~/blastPath/) { system("export BLASTDIR=$hit"); }
		if ($1=~/blastPath/) { system("export PATH=\$PATH:$hit"); }
		elsif ($1=~/blastDB/) { $blastParams[0]=$hit; }
		elsif ($1=~/threshold/) { $blastParams[1]=$hit; }
	#Now grab project primer parameters
	} elsif ($curLine=~m/projPrimer_(\w*) = (\d*)/) {
		my $hit=$2;
		if 	($1=~/lengthMin/) { 	$primer1Params[0]=$hit; }
		elsif 	($1=~/lengthMax/) { 	$primer1Params[1]=$hit; }	
		elsif 	($1=~/TmMin/) { 	$primer1Params[2]=$hit; }	
		elsif 	($1=~/TmMax/) { 	$primer1Params[3]=$hit; }	
		elsif 	($1=~/GCmin/) { 	$primer1Params[4]=$hit; }	
		elsif 	($1=~/GCmax/) { 	$primer1Params[5]=$hit; }	
		elsif 	($1=~/maxATCreps/) { 	$primer1Params[6]=$hit; }	
		elsif 	($1=~/maxGreps/) { 	$primer1Params[7]=$hit; }	
		elsif 	($1=~/minHbonds/) { 	$primer1Params[8]=$hit; }
		elsif 	($1=~/minPalLength/) { 	$primer1Params[9]=$hit; }
		elsif 	($1=~/maxPalLength/) { 	$primer1Params[10]=$hit; }
		elsif 	($1=~/gapLimit/) { 	$primer1Params[11]=$hit; }
		elsif 	($1=~/numMismatch/) { 	$primer1Params[12]=$hit; }
	#Now grab target primer parameters
	} elsif ($curLine=~m/targPrimer_(\w*) = (\d*)/) {
		my $hit=$2;
		if 	($1=~/lengthMin/) {	$primer2Params[0]=$hit; }
		elsif 	($1=~/lengthMax/) { 	$primer2Params[1]=$hit; }	
		elsif 	($1=~/TmMin/) { 	$primer2Params[2]=$hit; }	
		elsif 	($1=~/TmMax/) { 	$primer2Params[3]=$hit; }	
		elsif 	($1=~/GCmin/) { 	$primer2Params[4]=$hit; }	
		elsif 	($1=~/GCmax/) { 	$primer2Params[5]=$hit; }	
		elsif 	($1=~/maxATCreps/) { 	$primer2Params[6]=$hit; }	
		elsif 	($1=~/maxGreps/) { 	$primer2Params[7]=$hit; }	
		elsif 	($1=~/minHbonds/) { 	$primer2Params[8]=$hit; }
		elsif 	($1=~/minPalLength/) { 	$primer2Params[9]=$hit; }
		elsif 	($1=~/maxPalLength/) { 	$primer2Params[10]=$hit; }
		elsif 	($1=~/gapLimit/) { 	$primer2Params[11]=$hit; }
		elsif 	($1=~/numMismatch/) { 	$primer2Params[12]=$hit; }
	#Now grab target probes parameters
	} elsif ($curLine=~m/targProbe_(\w*) = (\d*)/) {
		my $hit=$2;
		if 	($1=~/lengthMin/) { 	$probeParams[0]=$hit; }
		elsif 	($1=~/lengthMax/) { 	$probeParams[1]=$hit; }	
		elsif 	($1=~/TmMin/) { 	$probeParams[2]=$hit; }	
		elsif 	($1=~/TmMax/) { 	$probeParams[3]=$hit; }	
		elsif 	($1=~/GCmin/) { 	$probeParams[4]=$hit; }	
		elsif 	($1=~/GCmax/) { 	$probeParams[5]=$hit; }	
		elsif 	($1=~/maxATCreps/) { 	$probeParams[6]=$hit; }	
		elsif 	($1=~/maxGreps/) { 	$probeParams[7]=$hit; }
		elsif 	($1=~/minHbonds/) { 	$probeParams[8]=$hit; }
		elsif 	($1=~/minPalLength/) { 	$probeParams[9]=$hit; }
		elsif 	($1=~/maxPalLength/) { 	$probeParams[10]=$hit; }
		elsif 	($1=~/gapLimit/) { 	$probeParams[11]=$hit; }
		elsif 	($1=~/numMismatch/) { 	$probeParams[12]=$hit; }	
	#Now grab spacer parameters
	} elsif ($curLine=~m/spacer_(\w*) = (\d*)/) {
		my $hit=$2;
		if 	($1=~/minHbonds/) { 	$spacerParams[0]=$hit; }
		elsif 	($1=~/minPalLength/) { 	$spacerParams[1]=$hit; }
		elsif 	($1=~/maxPalLength/) { 	$spacerParams[2]=$hit; }
		elsif 	($1=~/gapLimit/) { 	$spacerParams[3]=$hit; }
		elsif 	($1=~/numMismatch/) { 	$spacerParams[4]=$hit; }	
		elsif 	($1=~/primerGap/) { 	$spacerParams[5]=$hit; }
	} 
}

#Compile into single output w/ junk character spacers for later parsing
my @output=(@initParams,'x',@genParams,'x',@blastParams,'x',@primer1Params,
	'x',@primer2Params,'x',@probeParams,'x',@spacerParams,'x',
	@targetParams);
}
