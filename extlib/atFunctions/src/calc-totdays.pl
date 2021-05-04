#!/usr/local/bin/perl

@mdays0 = ( 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 );
@mdays1 = ( 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 );

$totday = 0;
print "{ ";
for ($i = 0; $i < 12; $i++) {
	print "${totday}, ";
	$totday += $mdays0[$i];
}
print "${totday} }\n";

$totday = 0;
print "{ ";
for ($i = 0; $i < 12; $i++) {
	print "${totday}, ";
	$totday += $mdays1[$i];
}
print "${totday} }\n";
