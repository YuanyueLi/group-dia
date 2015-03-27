#!/usr/bin/perl -w

open $fi,'<',$ARGV[0];
open $ft,'<',$ARGV[1];
open $fo,'>',$ARGV[2];

if($#ARGV!=2){
	die "iProphet_result  table_file Fixed_iProphet_result\n" 
}

while(<$ft>){
	chomp;
	@t=split "\t";
	$table{$t[0]}=$t[1];
}

while(<$fi>){
	if(/(^.+spectrum=".+?\.)\d+.\d+(.\d+" start_scan=")(\d+)(" end_scan=")\d+(".+$)/){
		$ori=$3;
		$scan=$table{$ori};
		print $fo $1.$scan.'.'.$scan.$2.$scan.$4.$scan.$5."\n";
	}else{
		print $fo $_;
	}
}
