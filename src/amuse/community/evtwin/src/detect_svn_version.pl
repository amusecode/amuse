#!/usr/bin/perl
use strict;
my $version;

if ($#ARGV>=0) {
   chdir $ARGV[0];
}

if (-d '.svn') {
   $version = `svnversion`;
}

if (defined $version && $version ne '') {
   $version = "svn r$version";
} else {
   $version = 'version unknown' 
}

chomp $version;
my $version_string = "TWIN/STARS ($version)";
my $n = length($version_string)+1;

open FILE, ">code/svn_version.f90";
print FILE "module svn_version
   character*$n :: svn_version_string = '$version_string'
end module
";
close FILE;
