#include <stdio.h>
#include <malloc.h>

#include "window.h"
#include "ntscDecode.h"

const unsigned int MAIN_WINDOW = 0;

const unsigned int MAIN_WIDTH = NTSC_SAMPLES_PER_LINE;
const unsigned int MAIN_HEIGHT = NTSC_LINES_PER_FRAME;

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

static void ntscEncodeTest(void)
{
	unsigned int x;
	unsigned int y;
	unsigned int i;

	const size_t frameSize = NTSC_LINES_PER_FRAME*NTSC_SAMPLES_PER_LINE;
	unsigned char* const output = malloc(frameSize);
	FILE* file = NULL;
	const char* const rgbFileName ="nes.rgb";
	const char* const ntscFileName ="jake.ntsc";
	size_t numBytesWritten = 0;

	const unsigned int NES_RGB_HEIGHT = 262;
	const unsigned int NES_RGB_WIDTH = 341;
	unsigned int* nesRGB = malloc(NES_RGB_HEIGHT*NES_RGB_WIDTH*sizeof(unsigned int));
	size_t numBytesRead = 0;

	file = fopen(rgbFileName, "rb");
	if (file == NULL)
	{
		fprintf(stderr, "ERROR can't open NES rgb data file '%s'\n", rgbFileName);
		return;
	}
	for (y = 0; y < NES_RGB_HEIGHT; y++)
	{
		for (x = 0; x < NES_RGB_WIDTH; x++)
		{
			unsigned int rgbValue;
			numBytesRead += fread(&rgbValue, 1, sizeof(unsigned int), file);
			nesRGB[y*NES_RGB_WIDTH+x] = rgbValue;
		}
	}
	printf("%d bytes read from '%s' width:%d height:%d\n", numBytesRead, rgbFileName, NES_RGB_WIDTH, NES_RGB_HEIGHT);
	fclose(file);

	ntscEncodeInit(output, frameSize);

	x = 0;
	y = 0;
	for (i = 0; i < frameSize; i++)
	{
		unsigned int pixelRGB = 0;
		unsigned int nesY = y/2 - 20;
		if (nesY < NES_RGB_HEIGHT)
		{
			unsigned int nesX = x/2 - 100;
			if (nesX < NES_RGB_WIDTH)
			{
				pixelRGB = nesRGB[nesY*NES_RGB_WIDTH+nesX];
			}
		}
		ntscEncodeAddSample(pixelRGB);
		x++;
		if (x >= NTSC_SAMPLES_PER_LINE)
		{
			x = 0;
			y += 2;
			if (y >= NTSC_LINES_PER_FRAME)
			{
				y -= NTSC_LINES_PER_FRAME;
			}
		}
	}
	file = fopen(ntscFileName, "wb");
	if (file == NULL)
	{
		fprintf(stderr, "ERROR can't open NTSC data file '%s'\n", ntscFileName);
		return;
	}
	numBytesWritten = fwrite((void*)output, 1, frameSize, file);
	if (numBytesWritten != frameSize)
	{
		fprintf(stderr, "ERROR writing frame data '%s' %d != %d\n", ntscFileName, numBytesWritten, frameSize);
		return;
	}
	fclose(file);
	free(output);
}

int main(int argc, char* argv[])
{
	int i;
	unsigned char* ntscDataPtr = NULL;
	unsigned int ntscDataSize = 0;
	unsigned int ntscFilenameIndex = 6;
	const char* ntscFilename = NULL;
	unsigned int* texture;

	char* ntscFilenames[] = {
		"smpte.ntsc",
		"lenna.ntsc",
		"indian_head.ntsc",
		"flightsim.ntsc",
		"penelope.ntsc",
		"nes_palette.ntsc",
		"jake.ntsc",
		NULL
	};

	for (i = 0; i < argc; i++)
	{
		printf("argv[%d] '%s'\n", i, argv[i]);
	}

	ntscEncodeTest();
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
