#!/usr/bin/perl -w
use Pg;


$conn = Pg::connectdb("dbname=expression_2");

$query = $conn->exec("select a.assembly, b.cstart, b.cend from fantom_assemblies a, fantom_matches b where a.assembly=b.assembly order by a.assembly");
while(@temp = $query->fetchrow){
    if(!defined($start{$temp[0]})){
	$start{$temp[0]} = $temp[1];
	$stop{$temp[0]} = $temp[2];
    }
    if($start{$temp[0]} > $temp[1]){ $start{$temp[0]} = $temp[1]; }
    if($stop{$temp[0]} < $temp[2]){ $stop{$temp[0]} = $temp[2]; }
}

for $id(sort {$a <=> $b} keys %start){
    print "$id\t$start{$id}\t$stop{$id}\n";
    $update = $conn->exec("update fantom_assemblies set begin=$start{$id} where assembly=$id");
    print "update Start : ", $update->cmdStatus , "\n";
    $update = $conn->exec("update fantom_assemblies set stop=$stop{$id} where assembly=$id");
    print "update Stop         : ", $update->cmdStatus , "\n";
}
