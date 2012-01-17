#include <stdio.h>
#include <malloc.h>

#include "window.h"

const unsigned int MAIN_WINDOW = 0;

const unsigned int MAIN_WIDTH = 910;
const unsigned int MAIN_HEIGHT = 262*2;

int loadNTSCData(const char* const fileName, unsigned char** pNtscDataPtr, unsigned int* pNtscDataSize)
{
	FILE* file;
	size_t size;
	size_t numBytesRead;

	*pNtscDataPtr = NULL;
	*pNtscDataSize = 0;

	file = fopen(fileName, "rb");
	if (file == NULL)
	{
		fprintf(stderr, "ERROR can't open NTSC data file '%s'\n", fileName);
		return 1;
	}
	fseek(file, 0, SEEK_END);
	size = (size_t)ftell(file);
	fseek(file, 0, SEEK_SET);

	*pNtscDataPtr = malloc(size);
	*pNtscDataSize = size;
	numBytesRead = fread((void*)*pNtscDataPtr, 1, size, file);
	if (numBytesRead != size)
	{
		fprintf(stderr, "ERROR can't reading NTSC data file '%s' numRead %d expected %d\n", fileName, numBytesRead, size);
		fclose(file);
		return 1;
	}
	fclose(file);

	printf("Loaded NTSC data file '%s' size %d bytes\n", fileName, size);
	return 0;
}

int main(int argc, char* argv[])
{
	int i;
	unsigned char* ntscDataPtr = NULL;
	unsigned int ntscDataSize = 0;

	if (loadNTSCData("2field-LLPPPP-110001.ntsc", &ntscDataPtr, &ntscDataSize))
	{
		fprintf(stderr, "ERROR loading NTSC data\n");
		return -1;
	}

	windowSetup(MAIN_WINDOW, MAIN_WIDTH, MAIN_HEIGHT);

	for (i = 0; i < argc; i++)
	{
		printf("argv[%d] '%s'\n", i, argv[i]);
	}
	while (1)
	{
		unsigned int x;
		unsigned int y;
		unsigned char* const texture = windowGetVideoMemoryBGRA(MAIN_WINDOW);
		unsigned int sample;
		sample = 0;
		i = 0;
		for (y = 0; y < MAIN_HEIGHT; y++)
		{
			for (x = 0; x < MAIN_WIDTH; x++)
			{
				if (1)
				{
					texture[i+0] = ntscDataPtr[sample];
					texture[i+1] = ntscDataPtr[sample];
					texture[i+2] = ntscDataPtr[sample];
					texture[i+3] = 0xFF;
					i += 4;
					sample++;
					if (sample >= ntscDataSize)
					{
						sample = 0;
					}
				}
				else
				{
					/* BGRA format */
					texture[i+0] = (unsigned char)((x * (y+10) + 100) & 0xFF);
					texture[i+1] = (unsigned char)(((x+10) * y + 50) & 0xFF);
					texture[i+2] = (unsigned char)(((x+5) * (y-10) + 20) & 0xFF);
					texture[i+3] = 0xFF;
					i += 4;
				}
			}
		}
		windowUpdate(MAIN_WINDOW, MAIN_WIDTH, MAIN_HEIGHT);
		windowMainLoop();

		/* 0x100 = ESC */
		if (windowCheckKey(0x100))
		{
			break;
		}
	}
	return -1;
}
