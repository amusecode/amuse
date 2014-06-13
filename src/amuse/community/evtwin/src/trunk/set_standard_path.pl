#!/usr/bin/perl
use strict;
my $version;

if ($#ARGV>=0) {
   chdir $ARGV[0];
}

my $datapath = "$ARGV[1]/stars";

$datapath =~ s/\/\//\//g;

my $n = length($datapath)+1;

if (-f "code/installation_path.f90") {
   open FILE, "<code/installation_path.f90";
   $_ = <FILE>;
   my $s = <FILE>;
   close FILE;

   chomp $s;
   $s =~ s/.*'(.*)'/$1/;

# Abort if nothing needs to be done, the file is up-to-date
   if ($s eq $datapath) { exit; }
}

open FILE, ">code/installation_path.f90";
print FILE "module install_location
   character(len=$n) :: twin_install_path = '$datapath'
end module
";
close FILE;
