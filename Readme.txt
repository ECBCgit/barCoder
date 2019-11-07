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
# Introduction
#########################################
# This document describes how to use the barcode generation code. The code will
# generate a list of sequences that can be directly synthesized and inserted
# into the target genomes via user-specified targeting vectors. 
#
# Running the code
#
# This software was tested on CentOS 7, and requires at least 64 GB of RAM in 
# order to support BLAST searches against the "nt" database. A modern CPU is 
# recommended, however the program is not currently optimized for multi-threaded
# processing. 
#
# Dependencies:
#
# * perl 5.8 or higher
# * Bioperl. Information on installing Bioperl can be found at
#   https://bioperl.org/INSTALL.html Specific Bioperl modules used in this program
#   include:
#	- StandAloneBlastPlus (bl2seq, blastn)
#	- EMBOSS
# * nt database from NCBI: ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*
#
# Test data is located in the Projects directory.	
#
# The purpose of the code is to generate a set of barcodes that minimize the chances of 
#		non-unique detection.
#
# Generating unique barcodes: 
# 	A project contains a set of genomes that need to have barcodes generated
# for them, called "target genomes". These barcodes meet a minimum set of 
# requirements for PCR-based detection, including standard PCR parameters
# (Tm, length, etc.) and high specificity. The specificity includes low homology
# to the target genomes, other primers, and a second set of genomes called
# "other genomes". This second set of genomes includes things that do not need
# to be barcoded, but are expected to be present in the environment. The primers
# must not match these other genomes in order to minimize false positives. The 
# primers are also blasted against the NCBI nt database to maximize the odds of
# uniqueness.
#	Three primers are generated. The first primer is the same across a 
# single project and is one of the two primers needed for qPCR detection using 
# either SYBR Green or TaqMan technologies. The second qPCR primer is unique for
# each target genome in the project. For each target genome a second unique 
# sequence is generated for use as an optional probe sequence, needed only for
# TaqMan-based qPCR detection.
#	A filler sequence completes the barcode. This sequence is randomly 
# generated with 2 constraints: (1) the GC content matches that of the target
# genome and (2) there are no major secondary structures to interfere with 
# the PCR. For each target genome, a second set of 3 primers are generated as
# a backup in case any of the sets fail for some unforseen reason.
#	At the time of execution, the user can choose the number of barcodes to 
# generate per target genome. This provides multiple alternative barcodes in the
# case that 1 or more do not function well.
#
#
#########################################
# User Guide
#########################################
# The first step is to create a new project directory. This directory needs to 
# 	be placed under the Barcode/Projects/ directory. There are directories
#	there called BAsterne, BTK and YPco92 as examples.
#
# This folder should contain the following (some are optional):
#	(1) Subdirectory containing the target genomes. These can be in genbank or
#		fasta format. 
#	(2) Text file containing the script parameters, called "parameters.txt".
#		These include things like target Tm and blast threshold to count 
#		as a "match." This file should be copied from the BAsterne, BTK or YPco92  
#		project and tweaked as appropriate. Each parameter is described 
#		in the "parameters.txt" file in the project directory.
#	(2) (optional) Subdirectory containing other genomes (genbank or fasta 
#		format).
#	(4) (optional) List of primers to avoid from previous projects. Fasta
#		file listing primers. File should be named prevPrimers* (i.e.
#		can have prevPrimers1.fasta, prevPrimersProject2.fas, etc.)
#
# Usage: 
# 	from the Barcoder/ directory, type:
#
#	perl barcoder.pl "project_name" nBarcodes
#
#	Where:
#	- project_name = name of the subdirectory of the current project, located
#		at Barcode/Projects/project_name
#	- nBarcodes = number of barcodes to generate per target genome (it is 
#		useful to generate multiple in case some fail). 
# Outputs:
#	- primerList.fa: list of primers generated for barcodes
#	- barcodeList.fa: list of barcodes generated
#	- log files: complete description of what the script did (see
#	 	parameters.txt for some options)
#
#########################################
