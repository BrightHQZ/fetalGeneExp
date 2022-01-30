#!/usr/bin/perl

use strict;
use Getopt::Long;

my ( $bedFile, $size, $center, $outF); 

GetOptions( 
   'bed:s'  => \$bedFile, 
   'bin:i' => \$size, 
   'center!' => \$center,
   'out:s' => \$outF,
); 

if (!(defined $bedFile) || !(-e $bedFile)) { print("The bed file is necessary in this program for --bed.\n"); exit(0); }
if (!(defined $size) || $size < 1) { print("The bin size is necessary in this program for --bin (size must over 1).\n"); exit(0); }
if (!(defined $outF)) { print("The output file is necessary in this program for --out.\n"); exit(0); }

open(rf, "<$bedFile") || die("File:$bedFile is not exist.\n");
open(wf, ">$outF") || die("File:$outF is not created.\n");

while(<rf>) {
    my @items = split("\t",$_);
    my @temp = @items;
    my @res;
    if ($items[2] - $items[1] < $size) { close(rf); close(wf); print("The bin size is larger than the range of splited!\n"); last; }
    
    if (defined $center) {
        my $offset = sprintf("%.0f",($size/2));
        for (my $i = $items[1] - $offset; $i <= $items[2]; $i=$i + $size) { $res[++$#res] = $i;}
        $res[0] = $res[0] + $offset;
        $temp[1] = $res[0]; $temp[2] = $res[0] + $offset; print wf join("\t",@temp);
        for (my $i = 1; $i < $#res; $i++) { $temp[1] = $res[$i]; $temp[2] = $res[$i] + $size; print wf join("\t",@temp); }
        $temp[1] = $res[$#res]; $temp[2] = $items[2]; print wf join("\t",@temp);
    } else {
        for (my $i = $items[1]; $i <= $items[2]; $i=$i + $size) { $res[++$#res] = $i;}
        foreach (@res) {
            if ($_ + $size <= $items[2]) { $temp[1] = $_; $temp[2] = $_ + $size; print wf join("\t",@temp); }
            else { $temp[1] = $_; $temp[2] = $items[2]; print wf join("\t",@temp); }
        }
    }
}
close(rf);
close(wf);