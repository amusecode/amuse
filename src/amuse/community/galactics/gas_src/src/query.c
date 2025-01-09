/*
 * This file contains four functions for prompting single value input, 
 * one for each of the data types char, int, float and double.
 * The functions are cquery, iquery, fquery and dquery.
 * Each function has two arguments: a prompt and a pointer variable of the
 * appropriate data type.  
 * The function will display the prompt (e.g. "Enter x") and a default
 * value.  The user can either enter a new value or accept the default
 * by hitting a return.
 * 
 * Example
 * main()
 * {
 *		char a[50];
 *		float x=5;
 *		
 *		strcpy(a,"hello there");
 *		cquery("Enter a phrase",a);
 *		fquery("Enter the value to be squared",&x);
 *		fprintf(stdout,"%s %g squared = %g\n",a,x*x);
 * }
 */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

void cquery(char *prompt, char *c)
{
	char *p, line[100], newprompt[200];
	fprintf(stderr, "%s [%s]: ", prompt, c);
	fgets(line, 80, stdin);
	if (*line != '\n')
	{
		for (p=line; *p != '\n' && *p != '#'; p++);
		if ( *p == '#' ) {
			p--;
			while ( *p == ' ' || *p == '\t' ) p--;
			p++;
		}
		*p = '\0';
		
		strcpy(c, line);
	}
}

void dquery(char *prompt, double *d)
{
	char *p, line[100];
	fprintf(stderr, "%s [%g]: ", prompt, *d);
	fgets(line, 80, stdin);
	if (*line != '\n')
	{
		for (p=line; *p !='\n' && *p !='#'; p++);
			*p = '\0';
		*d = atof(line);
	}
}

void fquery(char *prompt, float *f)
{
	char *p, line[100];
	fprintf(stderr, "%s [%g]: ", prompt, *f);
	fgets(line, 80, stdin);
	if (*line != '\n')
	{
		for (p=line; *p!='\n' && *p!='#'; p++);
			*p = '\0';
		*f = (float) atof(line);
	}
}

void iquery(char *prompt, int *i)
{
	char *p, line[100];
	fprintf(stderr, "%s [%d]: ", prompt, *i);
	fgets(line,80,stdin);
	if (*line != '\n')
	{
		for (p=line; *p !='\n' && *p !='#'; p++);
			*p = '\0';
		*i = atoi(line);
	}
}

void llquery(char *prompt, long long *i)
{
	char line[100];
	fprintf(stderr, "%s [%lld]: ", prompt, *i);
	fgets(line, 100, stdin);
	if (*line != '\n')
	{
		sscanf(line,"%lld", i);
	}
}

