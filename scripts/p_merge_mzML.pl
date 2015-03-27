#!/usr/bin/perl

use IO::File;
open $flist, '<', $ARGV[0];

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

# Find all mzML files.
@mzml_files = <$path*.mzML>;

@out;
for $i ( 0 .. ( $exp_size - 1 ) ) {
    $out[$i] = $name . '.spec.' . $i . '.mzML';
}
@fo = map { IO::File->new( $_, 'w' ) } @out;

#
# Go through the whole file, count spectrumList.
$spec_num = 0;
for $f (@mzml_files) {
    open $fi, '<', $f;
    while (<$fi>) {
        if (/<spectrum.+index="(.+)"\ defaultArrayLength/) {
            $scan  = $1;
            $decoy = $scan % 2;
            $scan  = ( $scan - $decoy ) / 2;
            $win   = $scan % $windows_num;
            $scan  = ( $scan - $win ) / $windows_num;
            $exp   = $scan % $exp_size;
            $spec_num[$exp]++;

        }

    }
}
#
# Merge all mzML files.
$out_head = 1;
$cur_scan = 1;
for $f (@mzml_files) {
    open $fi, '<', $f;
    $is_head = 1;
    $content = "";
    while (<$fi>) {
        $l = $_;
        if (/<spectrumList/) {
            if ($out_head) {
                for $i ( 0 .. ( $exp_size - 1 ) ) {
                    $fo[$i]->print(
                            '		<spectrumList count="'
                          . $spec_num[$i]
                          . '" defaultDataProcessingRef="dp_sp_0">
'
                    );
                }
            }
            $is_head  = 0;
            $out_head = 0;
            next;
        }

        if (/<\/spectrumList>/) {
            last;
        }
        if ( $is_head and $out_head ) {
            for $i ( 0 .. ( $exp_size - 1 ) ) {
                $fo[$i]->print($l);
            }
        }
        elsif ( !$is_head ) {
            if (
/(^.+spectrum\=)(\d+)(".+index=")(\d+)("\ defaultArrayLength.+$)/
              )
            {
				$origin_scan=$2;
                $scan  = $2;
                $decoy = $scan % 2;
                $scan  = ( $scan - $decoy ) / 2;
                $win   = $scan % $windows_num;
                $scan  = ( $scan - $win ) / $windows_num;
                $exp   = $scan % $exp_size;

                $content = $1 . $cur_scan . $3 . $cur_scan . $5 . "\n";
				$table{$cur_scan}=$origin_scan;
                $cur_scan++;
            }
            else {
                $content .= $_;
            }

            if (/<\/spectrum>/) {
                $fo[$exp]->print($content);
            }
        }
    }
}
for $fc (@fo) {
    $fc->print(
        '		</spectrumList>
	</run>
</mzML>'
    );
}
