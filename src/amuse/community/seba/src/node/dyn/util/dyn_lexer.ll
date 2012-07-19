/*
 *  Flex lexer to tokenize Starlab's dyn format.
 *  Copyright (C) 2003  StarCluster team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */


	#include <vector>
	#include "dyn.h"
	#include "dyn_parser.h"

	#if !HAVE_STPCPY
	  char *stpcpy(char *restrict, const char *restrict);
	#endif

	// we define this ourselves due to the peeking at the first character
	// in get_dyn() after which read(2) wouldn't work right
	#define YY_INPUT(buf, result, max_size) \
	  if ((result = fread(buf, 1, max_size, yyin)) < 0 ) \
	    YY_FATAL_ERROR("input in flex scanner failed");

	extern unsigned long lineno;
	extern struct { char* s; size_t n; } story_line;

%option noyywrap

%x LOG_TEXT LOG_ LOG__
%s VALUE


%%


Particle	return PARTICLE;
\(Log		yyless(0); BEGIN LOG__; return '(';
<LOG__>\(Log	BEGIN LOG_TEXT; return LOG;
<LOG_>\)Log	yyless(0); BEGIN LOG__; return ')';
<LOG__>\)Log	BEGIN INITIAL; return LOG;
Dynamics	return DYNAMICS;
Hydro		return HYDRO;
Star		return STAR;

<LOG_TEXT>.*\n							{
// log story: swallow everything up until ")Log"
  ++lineno;
  if (const char* const c = strstr(yytext, ")Log")) {
    // this is necessary because yyless is a macro that segfaults otherwise
    const size_t n = yyleng - strlen(c);
    yyless(n);
    --lineno;
    BEGIN LOG_;
  } else yytext[yyleng-1] = '\0';
  yylval.string = strdup(yytext);
  return LOG_STORY;
}

<VALUE>[[:alpha:]_(][[:alnum:]_(,)[:blank:]]*			{
// alphanumeric string with parentheses, underscores, commas, and whitespace
  yylval.string = strdup(yytext);
  BEGIN INITIAL;
  return STRING;
}

[[:alpha:]_][[:alnum:]_]*					{
// alphanumeric string with underscores
  yylval.string = strdup(yytext);
  BEGIN VALUE;
  return KEYWORD;
}

[-+]?([[:digit:]]*\.)?[[:digit:]]+([eE][-+]?[[:digit:]]+)?	{
// floating point number
  yylval.real = atof(yytext);
  story_line.s =
    static_cast<char*>(realloc(story_line.s, (story_line.n+=1+yyleng)+1));
  strcpy(stpcpy(story_line.s+story_line.n-yyleng-1, " "), yytext);
  return NUMBER;
}

\n		++lineno; BEGIN INITIAL;

[ \f\t\v]+	// skip whitespace

#.*		// skip the rest of line on comments

.		return yytext[0];
