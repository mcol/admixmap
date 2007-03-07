/* vim:set ft=lex: */

%{
/**
 * Lexer for the user data. It's job is to recognize strings which
 * appear in the input data.
 */

#include <stdio.h>
#include <string.h>

#include "UserData.h"
#include "UserDataIndividual.h"
#include "userdata-parser.h"

%}

%%

[ACGT]		yylval.text = strdup(yytext); return GENOTYPE_BASE;
:		return COLON;
rs[0-9]+	yylval.snp_name = strdup(yytext); return SNP_ID;
[0-9]+		yylval.number = atoi(yytext); return NUMBER;
\n		return LINE_END;
\t		return COLUMN_SEPARATOR;
\s		return COLUMN_SEPARATOR;
idno		return INDIVIDUAL_HEADER;
dncase		return OUTCOME_HEADER;
[a-z0-9_]+	yylval.text = strdup(yytext); return WORD;
%%
