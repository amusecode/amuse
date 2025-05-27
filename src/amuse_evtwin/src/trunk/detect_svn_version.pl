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

if (-f "code/svn_version.f90") {
   open FILE, "<code/svn_version.f90";
   $_ = <FILE>;
   my $s = <FILE>;
   close FILE;

   chomp $s;
   $s =~ s/.*'(.*)'/$1/;

# Abort if nothing needs to be done, the file is up-to-date
   if ($s eq $version_string) { exit; }
}

open FILE, ">code/svn_version.f90";
print FILE "module svn_version
   character(len=$n) :: svn_version_string = '$version_string'
end module
";
close FILE;
