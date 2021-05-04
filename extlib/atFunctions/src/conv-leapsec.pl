#!/usr/local/bin/perl

$leapsec_file = "/usr/local/astroe/com/calibration/leapsec/leapsec.fits";

open(fp, "fdump ${leapsec_file} STDOUT - - pagewidth=256 page=no prhead=no prdata=yes |");

while (<fp>) {
	($row, $date, $time, $mjd, $secs, $leapsecs) = (split);
	next if ( $row <= 0 );
	($yr,$mo,$dy) = (split(/-/, $date));
	($hr,$mn,$sc) = (split(/:/, $time));
	printf("\t{ { %04d, %02d, %02d, %02d, %02d, %02d, 0.0, 0.0 }, %.1f, %d },\n", $yr, $mo, $dy, $hr, $mn, $sc, $mjd, $leapsecs);
}
