#!/usr/bin/perl
use strict;

my $file;
my $moddep;

foreach $file (@ARGV) {
   my $basename = $file;
   my $objname;
   $basename =~ s/.*\///;

   $objname = $basename;
   $objname =~ s/\.(c|f)(90|95)?/\.o/;

   # Find module dependencies by looking for MODULE and USE statements
   my $modules;
   my @provides_modules;
   my @wants_modules;
   open SOURCE, "<$file";
   while (<SOURCE>) {
      chomp;
      if (m/^\s*use /i){
         my $module_name = $_;
         $module_name =~ s/\s*use\s*//i;
         $module_name =~ s/[^a-zA-Z0-9_]//g;
         $module_name =~ s/(\w+)/\L$1/g;
         unless (grep(/$module_name/,@wants_modules)) {
            push @wants_modules, $module_name;
         }
      }
      
      if (m/^\s*module /i){
         my $module_name = $_;
         $module_name =~ s/\s*module\s*//i;
         $module_name =~ s/[^a-zA-Z0-9_]//g;
         $module_name =~ s/(\w+)/\L$1/g;
         push @provides_modules, $module_name
      }
   }
   close SOURCE;

   map {
      # Don't depend on modules that are provided by this source file
      my $modname = $_;
      if (!grep /$modname/, @provides_modules) {
         $modules .= " \$(MODDIR)/$_.mod"
      }
   } (@wants_modules);
   map {
      $moddep  .= "\$(MODDIR)/$_.mod: \$(OBJDIR)/$objname\n";
   } (@provides_modules);
   undef @wants_modules;
   undef @provides_modules;
   
   print "\$(OBJDIR)/$objname: $file";
   print "$modules\n";
   if ($basename =~ m/\.c/) {
      print "\t\$(CC) \$(CFLAGS) \$(INCLUDE) -c -o \$@ \$<\n\n";
   } else {
      print "\t\$(FORT) \$(FFLAGS) \$(INCLUDE) -c -o \$@ \$<\n\n";
   }

   print "\$(OBJDIR)/pic_$objname: $file";
   print "$modules\n";
   if ($basename =~ m/\.c/) {
      print "\t\$(CC) \$(CFLAGS) -fPIC \$(INCLUDE) -c -o \$@ \$<\n\n";
   } else {
      #print "\t\$(FORT) \$(FFLAGS) \$(INCLUDE) -c -o \$@ \$<\n\n";
      print "\t\$(FORT) \$(FFLAGS) -fPIC \$(INCLUDE) -c -o \$@ \$<\n\n";
   }
}

print "$moddep";
