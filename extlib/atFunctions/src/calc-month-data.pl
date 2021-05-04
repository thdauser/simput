#!/usr/local/bin/perl

@mdays = (0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 );

print "static char month_norm[] = {\n";

for ($m = 1; $m <= 12; $m++) {
	$days = $mdays[$m];
	for ($d = 1; $d <= $days; $d++) {
		printf "%2d,", $m;
	}
	print "\t/* $days */\n";
}

print "0 };\n";
print "\n";
