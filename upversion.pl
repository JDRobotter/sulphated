#!/usr/bin/perl

$str = `cat ./version.h`;

($h) = ( $str =~ /#define VERSION_H (\d*)$/m );
($l) = ( $str =~ /#define VERSION_L (\d*)$/m );
($n) = ( $str =~ /#define VERSION_N (\d*)$/m );

$n = $n + 1;

print "#define VERSION_H $h\n";
print "#define VERSION_L $l\n";
print "#define VERSION_N $n\n";

