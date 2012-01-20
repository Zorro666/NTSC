#include <stdio.h>
#include <malloc.h>

#include "window.h"
#include "ntscDecode.h"

const unsigned int MAIN_WINDOW = 0;

const unsigned int MAIN_WIDTH = NTSC_SAMPLES_PER_LINE;
const unsigned int MAIN_HEIGHT = NTSC_LINES_PER_FIELD*NTSC_FIELDS_PER_IMAGE;

int loadData(const char* const fileName, unsigned char** pNtscDataPtr, unsigned int* pNtscDataSize)
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
	unsigned int ntscFilenameIndex = 0;
	const char* ntscFilename = NULL;
	unsigned int* texture;

	char* ntscFilenames[] = {
		"smpte.ntsc",
		"lenna.ntsc",
		"indian_head.ntsc",
		"flightsim.ntsc",
		"penelope.ntsc",
		"2field-LLPPPP-110001.ntsc",
		NULL
	};

	for (i = 0; i < argc; i++)
	{
		printf("argv[%d] '%s'\n", i, argv[i]);
	}

	ntscFilename = ntscFilenames[ntscFilenameIndex];
	if (loadData(ntscFilename, &ntscDataPtr, &ntscDataSize))
	{
		fprintf(stderr, "ERROR loading NTSC data '%s'\n", ntscFilename);
		return -1;
	}

	windowSetup(MAIN_WINDOW, MAIN_WIDTH, MAIN_HEIGHT);
	texture = windowGetVideoMemoryBGRA(MAIN_WINDOW);
	ntscDecodeInit(texture);

	while (1)
	{
		unsigned int sample;

		sample = 0;
		i = 0;
		for (i = 0; i < (int)(MAIN_HEIGHT*MAIN_WIDTH); i++)
		{
			const unsigned char sampleValue = ntscDataPtr[sample];
			ntscDecodeAddSample(sampleValue);
			sample++;
			if (sample >= ntscDataSize)
			{
				sample = 0;
			}
		}
		windowUpdate(MAIN_WINDOW, MAIN_WIDTH, MAIN_HEIGHT);
		windowMainLoop();
		ntscDecodeTick();

		/* 0x100 = ESC */
		if (windowCheckKey(0x100))
		{
			windowClearKey(0x100);
			break;
		}
		if (windowCheckKey('P'))
		{
			windowClearKey('P');
			ntscFilenameIndex++;
			ntscFilename = ntscFilenames[ntscFilenameIndex];
			if (ntscFilename == NULL)
			{
				ntscFilenameIndex = 0;
				ntscFilename = ntscFilenames[ntscFilenameIndex];
			}
			if (loadData(ntscFilename, &ntscDataPtr, &ntscDataSize))
			{
				fprintf(stderr, "ERROR loading NTSC data '%s'\n", ntscFilename);
				return -1;
			}
		}
	}
	return -1;
}
