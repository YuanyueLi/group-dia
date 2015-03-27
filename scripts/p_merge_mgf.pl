#!/usr/bin/perl -w
use IO::File;
open $flist,  '<', $ARGV[0];
open $ftable, '>', 'table';

#read list
$name = <$flist>;
chomp $name;
$windows = <$flist>;
chomp $windows;
$path = <$flist>;
<$flist>;
chomp $path;
$path .= '/';
$exp_size = 0;

while (<$flist>) {
    if (/mzML/) {
        $exp_size++;
    }
}

# read windows
open $fw, '<', $windows;
$windows_num = 0;
while (<$fw>) {
    chomp;
    if (/\d/) {
        $windows_num++;
    }
}

@file = <$path*.mgf>;
@out;
for $i ( 0 .. ( $exp_size - 1 ) ) {
    $out[$i] = $name . '.spec.' . $i . '.mgf';
}
@fo = map { IO::File->new( $_, 'w' ) } @out;
for $i ( 0 .. ( $exp_size - 1 ) ) {

    #open $fc, '>', $name . '.spec.' . $i . '.mgf';
    $fo[$i]->print( '
' );

    #    push @fo, $fc;
}

$cur_scan = 0;

for $file (@file) {
    $content = "";
    open $fi, '<', $file;
    while ( $l = <$fi> ) {
        if ( $l =~ /(^TITLE=.+\.)(\d+?)\.(\d+?)(\.\d+$)/ ) {
            $cur_scan++;
            $origin_scan      = $2;
            $l                = $1 . $cur_scan . '.' . $cur_scan . $4 . "\n";
            $table{$cur_scan} = $origin_scan;
        }
        if ( $l =~ /SCANS=(.+$)/ ) {
            $scan  = $1;
            $decoy = $scan % 2;
            $scan  = ( $scan - $decoy ) / 2;
            $win   = $scan % $windows_num;
            $scan  = ( $scan - $win ) / $windows_num;
            $exp   = $scan % $exp_size;

            $l = 'SCANS=' . $cur_scan . "\n";
        }

        if ( $l =~ /BEGIN\ IONS/ ) {
            $content = $l;
        }
        else {
            $content .= $l;
        }

        if ( $l =~ /END\ IONS/ ) {
            $fc = $fo[$exp];
            $fo[$exp]->print( $content . "\n" );
        }
    }
}

for $k ( keys %table ) {
    print $ftable $k . "\t" . $table{$k} . "\n";
}
