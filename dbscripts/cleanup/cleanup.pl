#!/usr/bin/perl -w
use Pg;

## this script generates an sql file which can be used to remove stuff from expression
## databases. Should remove..
## 
## 1. experiments not desired.
## 2. expression data not desired (as from above)
## 3. cel and other files for experiments not desired.
##

## Experiments Expression and Cel Files..
## read in a file containing the experiments we want to keep. 
## then contact database and remove those not included in there using the thingy..

if(@ARGV != 2){
    die "usage : cleanup.pl dbname retain_file\n";
}

($dbname, $keepFile) = @ARGV;

open(IN, $keepFile) || die "couldn't open $keepFile $!\n";
while(<IN>){
    chomp;
    @temp = split /\t/, $_;
    ## 0 index
    ## 1 experiment (double,,, treated like text anyway,, so)
    ## 2 expt_group
    ## 3 short_description
    ## 4 description
    
    ## include in a data structure so we can communicate what we will do to the user..
    $keepExpts{$temp[1]}{index} = $temp[0];
    $keepExpts{$temp[1]}{group} = $temp[2];
    $keepExpts{$temp[1]}{short} = $temp[3];
    $keepExpts{$temp[1]}{description} = $temp[4];
    $keeping{$temp[1]} = 1;   ## erase these when going through, so we can make sure the user isn't losing any
}

## next contact the database and check what experiments are there, then 
## work out which ones to discard. Then inform the user of this..
$conn = Pg::connectdb("dbname=$dbname");
$query = $conn->exec("select * from experiments");
while(@temp = $query->fetchrow){
    if(@temp){
	if(!defined($keepExpts{$temp[1]})){
	    $discardExpts{$temp[1]}{index} = $temp[0];
	    $discardExpts{$temp[1]}{group} = $temp[2];
	    $discardExpts{$temp[1]}{short} = $temp[3];
	    $discardExpts{$temp[1]}{description} = $temp[4];
	}else{
	    delete($keeping{$temp[1]});  ## variable name is not good. but allows us to check..stuff
	}
    }
}

## first let's check if the keep list contains anything that we couldn't find..
@missing = sort {$a <=> $b} keys %keeping;
print "Total number of experiments in keep list not found : ", $#keeping + 1, "\n";
for $expt(@missing){
    print "$keepExpts{$expt}{index}\t$keepExpts{$expt}{short}\t$expt\n";
}

## and then the number of experiments to delete..
@discardList = sort {$a <=> $b} keys %discardExpts;
print "Total number of experiments to delete : ", $#discardList + 1, "\n";
for $expt(@discardList){
    print "$discardExpts{$expt}{index}\t$discardExpts{$expt}{short}\t$expt\n";
}

## and then generate some sql script for this purpose.. 
## for each experiment that we need to delete we need to do three different things..

## 1. delete the data associated with the experiments (this takes a long time..)
## 2. delete cel files associated with the experiment .. 
## 3. delete the experiment row itself...
## these are the only tables that refer to the experiments as far as I can remember..

open(OUT, ">remove_experiments.sql") || die "couldn't open cleanup.sql $!\n";
print OUT "/*\nsql file created from $dbname use only with this database\n*/\n";
print OUT "begin;\n";
for $expt(@discardList){
    print OUT "delete from data where experiment=$expt;\n";
    print OUT "delete from cel_files where experiment=$expt;\n";
    print OUT "delete from experiments where experiment=$expt;\n";
}
print OUT "commit;\n";
