#!/usr/bin/perl -w
use Pg;

### Note this assumes that there is a rule on the user_annotation table that 
### removes any genes associated with a given comment. This isn't always the
### case, as I haven't yet decided whehter or not to make this a feature

if(@ARGV != 2){
    die "usage : cleanup_annotation.pl dbname retain_file\n";
}

($dbname, $keepFile) = @ARGV;

open(IN, $keepFile) || die "couldn't open $keepFile $!\n";
while(<IN>){
    @temp = split /\t/, $_;
    $keepIndex{$temp[0]} = 1;
}

$conn = Pg::connectdb("dbname=$dbname");
$query = $conn->exec("select * from user_annotation");
print "begin;\n";
while(@temp = $query->fetchrow){
    if(!defined($keepIndex{$temp[0]})){
	print "delete from user_annotation_genes where annotation_id = $temp[0];\n";
	print "delete from user_annotation where index=$temp[0];\n";
    }
}
print "commit;\n";
