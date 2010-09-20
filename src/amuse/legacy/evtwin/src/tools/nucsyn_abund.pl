#!/usr/bin/perl
use strict;

my %abund;
my %charge;

# Load abundance pattern from Anders & Grevesse
$abund{'XD'} = 4.801e-5;
$abund{'XHE3'} = 2.929e-5;
$abund{'XLI7'} = 9.353e-9;
$abund{'XBE7'} = 0;
$abund{'XB11'} = 4.725e-9;
$abund{'XC12'} = 3.520e-3;
$abund{'XC13'} = 3.650e-5;
$abund{'XC14'} = 0;
$abund{'XN14'} = 1.040e-3;
$abund{'XN15'} = 4.363e-6;
$abund{'XO16'} = 1.004e-2;
$abund{'XO17'} = 3.887e-6;
$abund{'XO18'} = 2.167e-5;
$abund{'XF19'} = 4.051e-7;
$abund{'XNE20'} = 1.84e-3;
$abund{'XNE21'} = 4.127e-6;
$abund{'XNE22'} = 1.302e-4;
$abund{'XNA22'} = 0;
$abund{'XNA23'} = 3.339e-5;
$abund{'XMG24'} = 5.148e-4;
$abund{'XMG25'} = 6.766e-5;
$abund{'XMG26'} = 7.760e-5;
$abund{'XAL26M'} = 0;
$abund{'XAL26G'} = 0;
$abund{'XAL27'} = 5.798e-5;
$abund{'XSI28'} = 6.530e-4;
$abund{'XSI29'} = 3.426e-5;
$abund{'XSI30'} = 2.352e-5;
$abund{'XP31'} = 8.155e-5;
$abund{'XS32'} = 3.958e-4;
$abund{'XS33'} = 3.222e-6;
$abund{'XS34'} = 1.866e-5;
$abund{'XFE56'} = 1.169e-3;
$abund{'XFE57'} = 2.855e-5;
$abund{'XFE58'} = 3.697e-6;
$abund{'XFE59'} = 0;
$abund{'XFE60'} = 0;
$abund{'XCO59'} = 3.358e-6;
$abund{'XNI58'} = 4.944e-5;
$abund{'XNI59'} = 0;
$abund{'XNI60'} = 1.958e-5;
$abund{'XNI61'} = 8.594e-7;

# Alternative: Grevesse & Noels (after Heger), supplemented by
# Anders&Grevesse
#      XD = 4.801d-5*(CZS/0.02)
#      XHE3 = 3.78d-5*(CZS/0.02)
#      XLI7 = 9.63d-9*(CZS/0.02)
#      XBE7 = 0d0
#      XB11 = 4.92d-9*(CZS/0.02)
#      XC12 = 3.62d-3*(CZS/0.02)
#      XC13 = 4.03d-5*(CZS/0.02)
#      XC14 = 0d0
#      XN14 = 1.07d-3*(CZS/0.02)
#      XN15 = 3.92d-6*(CZS/0.02)
#      XO16 = 1.04d-3*(CZS/0.02)
#      XO17 = 3.95d-6*(CZS/0.02)
#      XO18 = 2.08d-5*(CZS/0.02)
#      XF19 = 3.89d-7*(CZS/0.02)
#      XNE20 = 1.77d-3*(CZS/0.02)
#      XNE21 = 4.31d-6*(CZS/0.02)
#      XNE22 = 1.29d-4*(CZS/0.02)
#      XNA22 = 0d0
#      XNA23 = 3.60d-5*(CZS/0.02)
#      XMG24 = 5.60d-4*(CZS/0.02)
#      XMG25 = 7.09d-5*(CZS/0.02)
#      XMG26 = 7.81d-5*(CZS/0.02)
#      XAL26M = 0d0
#      XAL26G = 0d0
#      XAL27 = 6.24d-5*(CZS/0.02)
#      XSI28 = 7.08d-4*(CZS/0.02)
#      XSI29 = 3.58d-5*(CZS/0.02)
#      XSI30 = 2.37d-5*(CZS/0.02)
# The rest are from Anders & Grevesse
#      XP31 = 8.155d-5*(CZS/0.02)
#      XS32 = 3.958d-4*(CZS/0.02)
#      XS33 = 3.222d-6*(CZS/0.02)
#      XS34 = 1.866d-5*(CZS/0.02)
#      XFE56 = 1.169d-3*(CZS/0.02)
#      XFE57 = 2.855d-5*(CZS/0.02)
#      XFE58 = 3.697d-6*(CZS/0.02)
#      XFE59 = 0d0
#      XFE60 = 0d0
#      XCO59 = 3.358d-6*(CZS/0.02)
#      XNI58 = 4.944d-5*(CZS/0.02)
#      XNI59 = 0d0
#      XNI60 = 1.958d-5*(CZS/0.02)
#      XNI61 = 8.594d-7*(CZS/0.02)

# Isotope charge numbers
$charge{'XD'} = 1;
$charge{'XHE3'} = 2;
$charge{'XLI7'} = 3;
$charge{'XBE7'} = 4;
$charge{'XB11'} = 5;
$charge{'XC12'} = 6;
$charge{'XC13'} = 6;
$charge{'XC14'} = 6;
$charge{'XN14'} = 7;
$charge{'XN15'} = 7;
$charge{'XO16'} = 8;
$charge{'XO17'} = 8;
$charge{'XO18'} = 8;
$charge{'XF19'} = 9;
$charge{'XNE20'} = 10;
$charge{'XNE21'} = 10;
$charge{'XNE22'} = 10;
$charge{'XNA22'} = 11;
$charge{'XNA23'} = 11;
$charge{'XMG24'} = 12;
$charge{'XMG25'} = 12;
$charge{'XMG26'} = 12;
$charge{'XAL26M'} = 13;
$charge{'XAL26G'} = 13;
$charge{'XAL27'} = 13;
$charge{'XSI28'} = 14;
$charge{'XSI29'} = 14;
$charge{'XSI30'} = 14;
$charge{'XP31'} = 15;
$charge{'XS32'} = 16;
$charge{'XS33'} = 16;
$charge{'XS34'} = 16;
$charge{'XFE56'} = 26;
$charge{'XFE57'} = 26;
$charge{'XFE58'} = 26;
$charge{'XFE59'} = 26;
$charge{'XFE60'} = 26;
$charge{'XCO59'} = 27;
$charge{'XNI58'} = 28;
$charge{'XNI59'} = 28;
$charge{'XNI60'} = 28;
$charge{'XNI61'} = 28;


my $sum = 0.0;
map { $sum += $abund{$_} } (keys %abund);

print "&ABUND\n";
map {
   printf "   C%-6s = %8.3e\n", $_, $abund{$_}/$sum
} (sort {$charge{$a} <=> $charge{$b}} keys %abund);
print "/\n";

#map {
#  print "C$_, ";
#  printf "      DOUBLE PRECISION, SAVE :: C%-6s = %8.3e\n", $_, $abund{$_}/$sum
#} (sort {$charge{$a} <=> $charge{$b}} keys %abund);
