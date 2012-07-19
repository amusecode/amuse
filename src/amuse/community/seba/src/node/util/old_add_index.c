/*
 * Index Starlab data (.dyn) files for random access.
 * Copyright (C) 2005 Ernest N. Mamikonyan.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */


#include "config.h"
#include <getopt.h>
#include <string.h>
#include <stdbool.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <error.h>
#include <stdlib.h>
#include <inttypes.h>


size_t Snaps;
off_t *Idx, Curr_Char;
unsigned Depth;


static void usage(const int es) {
  fputs(
"Usage: add_index [OPTIONS] [FILE]\n"
"Index Starlab data (.dyn) files for random access.\n"
"Defaults are shown in parentheses.\n\n"
" -A, --append		append the index to file,\n"
"			if input is a stream, output just the index;\n"
"			this is the fastest mode of operation\n"
" -h, --help		display this help message and exit\n"
" -o, --output=FILE	place output in FILE (stdout)\n"
" -V, --version		display version and copyright information and exit\n",
       es == EXIT_SUCCESS ? stdout : stderr);
  exit(es);
}

static inline char proc_char(const char c) {
  if (c == '(') {
    if (!Depth++) {
      if (!(Idx = realloc(Idx, ++Snaps*sizeof(*Idx))))
	error(EXIT_FAILURE, errno,
	      "can't allocate sufficient memory: realloc(3)");
      Idx[Snaps-1] = Curr_Char;
    }
  } else if (c == ')') --Depth;
  return c;
}


int main(int argc, char *argv[]) {

  bool append = false;
  char c;
  size_t i;

  while (true) {
    int option, index;
    const struct option long_options[] = {
      { "append", no_argument, 0, 'A' },
      { "help", no_argument, 0, 'h' },
      { "output", required_argument, 0, 'o' },
      { "version", no_argument, 0, 'V' },
      { 0, 0, 0, 0 }
    };
    option = getopt_long(argc, argv, "Aho:V", long_options, &index);
    if (option == -1) break;
    switch (option) {
    case 'A': append = true; break;
    case 'h': usage(EXIT_SUCCESS);
    case 'o':
      if (strcmp(optarg, "-") && !(stdout = freopen(optarg, "w", stdout)))
	error(EXIT_FAILURE, errno, "can't open `%s' for writing: freopen(3)",
	      optarg);
      break;
    case 'V':
      puts(
"add_index (" PACKAGE_STRING ") 0.3\n"
"Report bugs to <" PACKAGE_BUGREPORT ">\n"
"Written by Ernest N. Mamikonyan.\n\n"
"Copyright © 2005 Starlab development group.\n"
"This is free software; see the source for copying conditions.  There is NO\n"
"warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.");
      exit(EXIT_SUCCESS);
    default: usage(EXIT_FAILURE);
    }
  }

  if (optind < argc) {
    if (!(stdin = freopen(argv[optind], "r", stdin)))
    error(EXIT_FAILURE, errno, "can't open `%s' for reading: freopen(3)",
	  argv[optind]);
    if (append && !(stdout = freopen(argv[optind], "a", stdout)))
      error(EXIT_FAILURE, errno, "can't open `%s' for writing: freopen(3)",
	    argv[optind]);
  }

#ifdef HAVE_MMAP
  struct stat statbuf;
  char *data;

  if (fstat(STDIN_FILENO, &statbuf) == -1)
    error(EXIT_FAILURE, errno, "can't stat `%s': fstat(2)", argv[optind]);

  if ((data = mmap(0, statbuf.st_size, PROT_READ, MAP_SHARED, STDIN_FILENO, 0))
      == MAP_FAILED)
    if (errno == ENOMEM || errno == ENODEV) goto read_with_stdio;
    else
      error(EXIT_FAILURE, errno, "can't map `%s' to memory: mmap(2)",
	    argv[optind]);

  if (append)
    while (Curr_Char < statbuf.st_size) proc_char(data[Curr_Char]), ++Curr_Char;
  else
    while (Curr_Char < statbuf.st_size)
      putchar_unlocked(proc_char(data[Curr_Char])), ++Curr_Char;

  goto print_idx;
#endif

 read_with_stdio:
  if (append)
    while ((c = getchar_unlocked()) != EOF) proc_char(c), ++Curr_Char;
  else
    while ((c = getchar_unlocked()) != EOF)
      putchar_unlocked(proc_char(c)), ++Curr_Char;

 print_idx:
  if (Depth)
    error(EXIT_FAILURE, 0, "premature EOF encountered or data is corrupted");

  fputs("index:", stdout);
  for (i = 0; i < Snaps; printf(" %" PRId64, Idx[i++]));
  putchar('\n');

}
