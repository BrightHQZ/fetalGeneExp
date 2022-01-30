#!/usr/bin/perl

use strict;
use Getopt::Long;

my ( $bedFile, $outD); 

GetOptions( 
   'bed:s'  => \$bedFile, 
   'out:s' => \$outD,
); 

if (!(defined $bedFile) || !(-e $bedFile)) { print("The bed file is necessary in this program for --bed.\n"); exit(0); }
if (!(defined $outD)) { print("The output dir is necessary or can't created in this program for --out $outD.\n"); exit(0); }

if (!(-d $outD)) {
    mkdir($outD) || die("Creating dir: $outD fail!");
} 

open(rf, "<$bedFile") || die("File:$bedFile is not exist.\n");
my %index;

while(<rf>) {
    $_ =~ s/\n//;
    my @items = split("\t",$_);
    #my @gene = split("_", $items[3]);
    if (!exists($index{$items[3]})) {
        $index{$items[3]} = 1;
        open(wf, ">".$outD."/".$items[3].".bed") || die("Can't create file:$outD"."\\".$items[3].".bed");
        print wf $items[0]."\t".$items[1]."\t".$items[2]."\t".$items[5]."\t".$items[4]."\t".$items[3]."\n";
    } else {
        print wf $items[0]."\t".$items[1]."\t".$items[2]."\t".$items[5]."\t".$items[4]."\t".$items[3]."\n";
    }
}
close(rf);
close(wf);