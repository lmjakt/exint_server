#!/usr/bin/perl -w
use Pg;

## sort out the thingy table..

$conn=Pg::connectdb("dbname=expression_2");
$query = $conn->exec("select * from p_sets");
while(@temp = $query->fetchrow){
    if(@temp){
	$afid{$temp[0]} = $temp[3];
	$id{$temp[0]} = $temp[4];
	$uid{$temp[0]}{$temp[6]} = 1;
	$uid{$temp[0]}{$temp[9]} = 1;
	$description{$temp[0]} = $temp[7];
    }
}

$query = $conn->exec("select af_id, description from tigr_annotation");
while(@temp = $query->fetchrow){
    if(@temp){
	$tigr{$temp[0]} = $temp[1];
    }
}

$query = $conn->exec("select * from uni_data");
while(@temp = $query->fetchrow){
    if(@temp){
	$uni{$temp[0]}{title} = $temp[1];
	$uni{$temp[0]}{gene} = $temp[2];
    }
}

## and put it all together..
for $index(sort {$a <=> $b} keys %afid){
    for $uniid(keys %{$uid{$index}}){
	## print out one line for each of these.. 
	if(defined($uni{$uniid}{gene})){
	    print "$index\t$afid{$index}\t$id{$index}\t$uniid\t$uni{$uniid}{gene}\t$uni{$uniid}{title}\t$description{$index}\t";
	    if(defined($tigr{$afid{$index}})){
		print $tigr{$afid{$index}};
	    }else{
		print "";
	    }
	    print "\n";
	}
    }
}
    
	     
