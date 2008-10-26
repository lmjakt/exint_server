#!/usr/bin/perl -w
use Pg;

## this script generats sql commands which can be piped into an expression database
## allowing the removal of data from in_situ_probes, and tables linked to these 
## i.e. images, annotation, image annotation, and others.

if(@ARGV != 2){
    die "useage : cleanup_ish_probes.pl dbname retain_file\n";
}

($dbname, $keepFile) = @ARGV;

### first read in the data in the keep file. This file is just a dump of the 
### the in_situ_probes table, with the rows to keep present, and other rows deleted.

open(IN, $keepFile) || die "couldn't open $keepFile $!\n";
while(<IN>){
    @temp = split /\t/, $_;  ## 0 is the index, which is the only thing we need really..
    $keepIds{$temp[0]} = 1;
}

### then work out which ones to delete..
$conn = Pg::connectdb("dbname=$dbname");
$query = $conn->exec("select probe, probe_id from in_situ_probes");
$ntuples = $query->ntuples;
print "Got a total of $ntuples tuples\n";
open(OUT, ">remove_ish_probes.sql") || die "couldn't open remove_ish_probes.sql\n";
print OUT "/*\nRemoves in situ probes. Most of the work is actually implemented in a set of rules in the database, which cause images and annotation associated with a probe to be removed as well\nwell it should be but there seems to be some bug with pg_dump, pg_restore and so on.. so most will be manua.\n*/\n";
print OUT "begin;\n";
while(@temp = $query->fetchrow){
    if(!defined($keepIds{$temp[0]})){
	### work out all the things we need to delete..
	## we need to delete from ..
	## 1. ish_probe_blast  
	## 2. ish_probe_assemblies
	##    a. ish_probe_matches
	## 3. ish_probe_classification
	## 4. ish_probe_num_annotation
	## 5. ish_probe_text_annotation
	## 6. ish_images
	##      ---- Note requires separate unlinking of images -- have to work out how to do this.
	##      ---- This is now handled by a rule defined as .. 
	##    dropimage as on delete to ish_images do (select lo_unlink(old.image); select lo_unlink(old.thumbnail);)
        ##    a. ish_annotation
	##
	## This is not too bad, but it doesn't touch the experiments or protocol
	## tables. However, it is not exactly clear how those should be handled.
	## 
	## Note that it is possible to set up rules for everything, such that a simple delete from ish_probes
	## would do everything. ----- but.... hmm, this seems like it might be a waste of time.. 
	## OK, this has now been done, and it should be pretty quick to remove the 
	## the comments..
	
	## the blast..
#	print "index : $temp[0]\t$temp[1]\n";
	print OUT "delete from ish_probe_blast where probe=$temp[0];\n";
	## then obtain the assembly ids from the the assembly table
	$query2 = $conn->exec("select assembly from ish_probe_assemblies where probe=$temp[0]");
	while(@temp2 = $query2->fetchrow){
	    print OUT "delete from ish_probe_matches where assembly=$temp2[0];\n";
	}
	print OUT "delete from ish_probe_assemblies where probe=$temp[0];\n";
	print OUT "delete from ish_probe_text_annotation where probe_id=$temp[0];\n";
	print OUT "delete from ish_probe_num_annotation where probe_id=$temp[0];\n";
	print OUT "delete from ish_probe_classification where probe_id=$temp[0];\n";
	### and then from images. We need to remove image annotation first..
	##  it seems I can still dump stuff with this image on it. Strange things are happening at the moment though.
	print OUT "delete from ish_images where probe=$temp[0];\n";
	print OUT "delete from in_situ_probes where probe=$temp[0];\n";
    }
}
print OUT "commit;\n";

