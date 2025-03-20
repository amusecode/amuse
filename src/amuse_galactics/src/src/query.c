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
 *      fprintf(stdout,"%s %g squared = %g\n",a,x*x);
 * }
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void cquery(char *prompt, char *c)
{
	char *p, line[100];
	fprintf(stderr, "%s [%s]: ", prompt, c);
	fgets(line, 100, stdin);
	if(*line != '\n')
	{
		/* overwrite what the user just typed with zeros and then copy?!! */
		for(p=line; *p != '\n'; p++);
			*p = '\0';
		strcpy(c, line);
	}
}

void dquery(char *prompt, double *d)
{
	char line[100];
	fprintf(stderr, "%s [%g]: ", prompt, *d);
	fgets(line, 100, stdin);
	if(*line != '\n')
	{
		sscanf(line, "%lf", d);
	}
}

void fquery(char *prompt, float *f)
{
	char line[100];
	fprintf(stderr, "%s [%g]: ", prompt, *f);
	fgets(line, 100, stdin);
	if(*line != '\n')
	{
		sscanf(line, "%e", f);
	}
}

void iquery(char *prompt, int *i)
{
	char line[100];
	fprintf(stderr, "%s [%d]: ", prompt, *i);
	fgets(line, 100, stdin);
	if(*line != '\n')
	{
		sscanf(line, "%d", i);
	}
}

void llquery(char *prompt, long long *i)
{
	char line[100];
	fprintf(stderr, "%s [%lld]: ", prompt, *i);
	fgets(line, 100, stdin);
	if(*line != '\n')
	{
		sscanf(line, "%lld", i);
	}
}

