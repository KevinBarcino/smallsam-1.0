// line 663

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{
	FILE *pFile1 = NULL;
	FILE *pFile2 = NULL;
	FILE *pFile3 = NULL;
	int c = 0;
	int count = 0;
	int line = 0;
	int n = 0;
	int i = 0, j = 0;
        int x = 0, y = 0;
	char *data = NULL;
	char buffer[10];
	char filename[2048] = {0};
	char fileSorted[2048] = {0};
	char fileData[2048] = {0};
	char filenameBuf[2048] = {0};
	
	if (argc == 1)
	{
		printf("Missing file name: e.g. ex1 filename\n");
		return 0;
	}
	else
	{	
		char *pos = NULL;
		strcpy(filenameBuf, argv[1]);
		strcpy(filename, filenameBuf);
		pos = strrchr(filenameBuf, '.');
		if (pos != NULL)
		{
			*pos = '\0';
		}	
		strcpy(fileSorted, filenameBuf);
		strcat(fileSorted, "-sorted.txt");
		strcpy(fileData, filenameBuf);
		strcat(fileData, "-data.csv");		

		printf("%s\n", filename);
		printf("%s\n", fileSorted);
		printf("%s\n", fileData);
		
	}

	// Find the number of tabs
	pFile1 = fopen(filename, "r");
	if (pFile1 == NULL)
	{	
		perror("Error opening file(s)");
		return 1;	
	}
	else
	{
		c = fgetc(pFile1);
		while (c != EOF)
		{	
			if (c == '\t' || c == '\x0A')
			{
				++count;
			}
			c = fgetc(pFile1);
		}
		fclose(pFile1);
	}

	// Allocate memory
	data = (char *)malloc(count*100);
	if (data == NULL)
	{
		perror("Error allocating memory.");
		return 1;
	}

	// put each element in an array
	for (i = 0; i < count; i++)
		*(data + i * 100) = '\0';
	pFile1 = fopen(filename, "r");
	pFile2 = fopen(fileSorted, "w");
	pFile3 = fopen(fileData, "w");
	//pFile1 = fopen("toy.sam", "r");
	//pFile2 = fopen("toy-sorted.txt", "w");
	//pFile3 = fopen("toy-data.csv", "w");
	char dataBuffer[1000] = {0};
	char data3[3][1000];
	char buf[100];
	char *pos = NULL;
	int seqCount = 0;
	int beginPos = 0;
	int length = 0;
	int ignoreThisLine = 0;

	if (pFile1 == NULL || pFile2 == NULL || pFile3 == NULL)
	{
		perror("Error opening file(s)");
	}
	else
	{
		sprintf(dataBuffer, "Begin\tEnd(M)\tLen\tEnd\tSequence\n");
		fputs(dataBuffer, pFile3);
		x = 0;
		y = 0;
		c = fgetc(pFile1);
		while (c != EOF)
		{	
			if (c == '@')
				ignoreThisLine = 1;
			else if ((c == '\t' || c == '\x0A'))
			{
				*(data + x * 100 + y) = '\0';

				sprintf(buffer, "%05d ", line + 1);
				// write sorted data to file
				fputs(buffer, pFile2);
				fputs((data + x *100), pFile2);
				fputs("\n", pFile2);

				if (c == '\t')
				{

					// write data to file
					if (seqCount == 2 || seqCount == 3 || seqCount == 5)
					{ 
						if (seqCount == 2)
							strcat(data3[0],data + x * 100);
						else if (seqCount == 3)
							strcat(data3[1],data + x * 100);
						else if (seqCount == 5)
							strcat(data3[2],data + x * 100);
						if (!ignoreThisLine && seqCount == 5)
						{
							strcpy(buf, data3[2]);
							if ((pos = strchr(data3[2], 'M')) != NULL)
							{
								*(buf + (pos - data3[2])) = '\0';
							}
							else
							{
								strcpy(buf, "0");
							}
							beginPos = atoi(data3[1]);
							length = atoi(buf);
							sprintf(dataBuffer, "%s\t%s\t%d\t%d\t%s\n", data3[1], data3[2], length, beginPos+length, data3[0]);
							fputs(dataBuffer, pFile3);
						}
					}
					line++;
					seqCount++;
				}
				else if (c == '\x0A')
				{
					fputs("**************************************************\n", pFile2);
					//fputs("\n", pFile3);
					seqCount = 0;
					*dataBuffer = '\0';
					*data3[0] = '\0';
					*data3[1] = '\0';
					*data3[2] = '\0';
					ignoreThisLine = 0;
					line++;
				}
				x++;
				y = 0;
			}
			else
			{
				*(data + x * 100 + y) = c;
                                y++;
			}
			c = fgetc(pFile1);
			
		}
		fclose(pFile1);
		fclose(pFile2);
		fclose(pFile3);

		for (i = 0; i < count; i++)
			if (*(data + i*100) != '\0')
				printf("%5d %s***\n", i + 1, (data + i * 100));

		printf("The file contains %d tab characters.\n", count);
	}

	free(data);
	return 0;
}



