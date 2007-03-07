/**
 * itoa() from http://www.freebookzone.com/others/itoa.h
 */
char *itoa(int n)
{
	int i = 0, j;
	char *s;
	char *u;

	s = (char *) malloc(17);
	u = (char *) malloc(17);

	do {
		s[i++] = (char) (n % 10 + 48);
		n -= n % 10;
	}
	while ((n /= 10) > 0);
	for (j = 0; j < i; j++)
		u[i - 1 - j] = s[j];

	u[j] = '\0';
	return u;
}

