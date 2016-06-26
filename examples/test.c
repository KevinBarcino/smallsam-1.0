#include <stdio.h>
int main() {
	FILE *pFile1, *pFile2;
	int c;
	int count = 0;
	int n = 0;
	char data[2048][100];
	int i, j;
	for (i = 0; i < 2048; i++)
		for(j = 0; j < 100; j++)
			data[i][j] = '0';
	pFile1 = fopen("ex1.sam", "r");
	pFile2 = fopen("ex1b.sam", "w");
	if (pFile1 == NULL || pFile2 == NULL)
	{
		perror("Error opening file(s)");
	}
	else
	{
		c = fgetc(pFile1);
		while (count < 5000 && c != EOF)
		{
			fputc(c, pFile2);
			c = fgetc(pFile1);
			++count;
		}
		fclose(pFile1);
		fclose(pFile2);
		printf("The file contains %d characters ($).\n", count);
	}
	return 0;
}
