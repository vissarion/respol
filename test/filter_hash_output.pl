#!/usr/bin/perl

# filter the hash output and show only how many homogeneous determinants
# were computed (and hashed) and how much times the determinant was called

while($line=<STDIN>){
        if ($line =~ m/non-hom d/ ){
                chomp $line;
                $line =~ s/.*\ (\d+)\ out\ of\ (\d+)/$1 $2/;
                print $line;
        }
}
