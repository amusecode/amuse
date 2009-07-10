#!/usr/bin/perl
# Configure sript for the TWIN stellar evolution code and library
# Pulled into MUSE for use with the SWIG wrapper around SSE
use strict;

# First: detect the name of the system
# Needed because the linker options are different in Linux than in OS X
print "Determining system name... ";
my $uname;
open(UNAME, "uname|") || die "Error: cannot run uname: $!";
$uname = <UNAME>;
chomp $uname;
close UNAME;
print "$uname\n";

my %platform_options;

# Use either the C or FORTRAN compiler in the linking stage; this is so it
# will pull in the correct dependencies and search the right directories
if ($uname =~ /Darwin/i) {
   $platform_options{'LDFLAGS'} = "-Wl,-single_module -dynamiclib -Wl,-flat_namespace -Wl,-undefined -Wl,suppress";
   $platform_options{'LD'} = "\$(CC)";
   $platform_options{'LIBS'} = "\$(FLIB) -lm";
   $platform_options{'SONAME'} = "dylib";
} else {
   $platform_options{'LDFLAGS'} = "-shared";
   $platform_options{'LD'} = "\$(FORT)";
   $platform_options{'LIBS'} = "\$(FLIB)";
   $platform_options{'SONAME'} = "so";
}

# Detect Python and Python version
print "Looking for Python... ";
my $python = `which python`;
chomp $python;
if (!$python) {
   print "not installed\n";
} else {
   print "$python\n";
   $platform_options{'PYTHON'} = "$python";

   print "Determining Python version... ";
   my $pyver = `$python -c "import sys; print (\\\"%d.%d\\\" % (sys.version_info[0], sys.version_info[1]))"`;
   chomp $pyver;
   print "$pyver\n";

   print "Looking for Python.h... ";
# First: try getting the information directly from Python
   my $python_root = `$python -c "import sys; print sys.prefix"`;
   chomp $python_root;
   my $python_path = "$python_root/include/python$pyver";
   if (-f "$python_path/Python.h") {
      print "$python_path/Python.h\n";
      $platform_options{'PYTHON_INCLUDE_DIR'} = $python_path;
   } else {
# Doesn't seem to be the correct path, try locate
      my $python_path = `locate python$pyver/Python.h`;
      if ($python_path) {
         split "\n", $python_path;
         $python_path = $_[0];
         print "$python_path\n";
         $python_path =~ s/Python.h$//;
         $platform_options{'PYTHON_INCLUDE_DIR'} = $python_path;
      } else {
         print "not found\n";
      }
   }

# Determine platform type of Python executable (the library we build will
# have to amtch this)
   print "Detecting architecture type...";
   my $platform_string = `file $python`;
   if ($platform_string =~ /64-bit/) { # 64-bit Python
      print " 64 bit\n";
      $platform_options{'ARCHITECTURE_FLAG'} = '-m64';
   } elsif ($platform_string =~ /32-bit/ || $platform_string =~ /i386/) { # 32-bit Python
      print " 32 bit\n";
      $platform_options{'ARCHITECTURE_FLAG'} = '-m32';
   } else { # No clue, assume default setting will do
      print " unknown\n";
      $platform_options{'ARCHITECTURE_FLAG'} = '';
   }
}

print "\nDetected these options:\n";
open MAKEFILE, ">makefile.sys";
foreach (keys %platform_options) {
   print "$_ = $platform_options{$_}\n";
   print MAKEFILE "$_ = $platform_options{$_}\n";
}
close MAKEFILE;
