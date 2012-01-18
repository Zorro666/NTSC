#include <stdio.h>
#include <math.h>
#include <malloc.h>

#include "window.h"

const unsigned int MAIN_WINDOW = 0;

const unsigned int MAIN_WIDTH = 910;
const unsigned int MAIN_HEIGHT = 262*2;

/*
1. [gain, level] = normalize(input)
2. input = input ∗ gain + level
4. Y = LPF2(video)
5. Sin = gen sin ( color burst(input))
6. I = LPF3( (video − Y ) . ∗ Sin )
7. Q = LPF3( (video − Y ) . ∗ Sin(2 :) )
8. [R G B]ˆT = clip(M−1 ∗ [Y I Q]ˆT
*/

/*
r  = rez amount, from sqrt(2) to ~ 0.1
f  = cutoff frequency
(from ~0 Hz to SampleRate/2 - though many synths seem to filter only  up to SampleRate/4)

The filter algo:
out(n) = a1 * in + a2 * in(n-1) + a3 * in(n-2) - b1*out(n-1) - b2*out(n-2)

Lowpass:
      c = 1.0 / tan(pi * f / sample_rate);

      a1 = 1.0 / ( 1.0 + r * c + c * c);
      a2 = 2* a1;
      a3 = a1;
      b1 = 2.0 * ( 1.0 - c*c) * a1;
      b2 = ( 1.0 - r * c + c * c) * a1;

Hipass:
      c = tan(pi * f / sample_rate);

      a1 = 1.0 / ( 1.0 + r * c + c * c);
      a2 = -2*a1;
      a3 = a1;
      b1 = 2.0 * ( c*c - 1.0) * a1;
      b2 = ( 1.0 - r * c + c * c) * a1;
*/

/* pSamples = input[n-2] */
/* pOutput = output[n-2] */
void computeLowPassCoeffs(float a[3], float b[2], const float freq, const float sample_rate)
{
	const float r = 1.0f;
	const float c = 1.0f / tanf((float)(M_PI * freq/sample_rate));

	a[0] = 1.0f / (1.0f + r*c + c*c);
	a[1] = 2.0f * a[0];
	a[2] = a[0];
	b[0] = 2.0f * (1.0f - c*c) * a[0];
	b[1] = (1.0f - r*c + c*c) * a[0];
}

void lowPass( const float a[3], const float b[2], const float* const pSamples, float* const pOutput)
{
	/* out(n) = a1 * in + a2 * in(n-1) + a3 * in(n-2) - b1*out(n-1) - b2*out(n-2) */
	pOutput[2] = a[0]*pSamples[2] + a[1]*pSamples[1] + a[2]*pSamples[0] - b[0]*pOutput[1] - b[1]*pOutput[0];
}

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
	int displayMode = 0;
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
		unsigned int Y = 0;
		sample = 0;
		i = 0;
		for (y = 0; y < MAIN_HEIGHT; y++)
		{
			for (x = 0; x < MAIN_WIDTH; x++)
			{
				const unsigned char sampleValue = ntscDataPtr[sample];
				unsigned char red = 0;
				unsigned char green = 0;
				unsigned char blue = 0;
				unsigned char alpha = 0xFF;
				if (displayMode == 0)
				{
					red = sampleValue;
					green = sampleValue;
					blue = sampleValue;
				}
				else if (displayMode == 1)
				{
				 	/* Y = sampleValue * (B-A) + y * B * A */
					const unsigned int SCALE = 1000;
					const unsigned int A = 990;
					const unsigned int B = SCALE - A;
				 	Y = (unsigned int)(sampleValue * B + Y * A);
					Y /= SCALE;
					red = (unsigned char)Y;
					green = (unsigned char)Y;
					blue = (unsigned char)Y;
				}
				else if (displayMode == 99)
				{
					red = (unsigned char)(((x+5) * (y-10) + 20) & 0xFF);
					green = (unsigned char)(((x+10) * y + 50) & 0xFF);
					blue = (unsigned char)((x * (y+10) + 100) & 0xFF);
				}
				/* BGRA format */
				texture[i+0] = blue;
				texture[i+1] = green;
				texture[i+2] = red;
				texture[i+3] = alpha;
				i += 4;
				sample++;
				if (sample >= ntscDataSize)
				{
					sample = 0;
				}
			}
		}
		windowUpdate(MAIN_WINDOW, MAIN_WIDTH, MAIN_HEIGHT);
		windowMainLoop();

		/* 0x100 = ESC */
		if (windowCheckKey(0x100))
		{
			windowClearKey(0x100);
			break;
		}
		/* 0x44 = d */
		if (windowCheckKey(0x44))
		{
			displayMode++;
			if (displayMode > 5)
			{
				displayMode = 0;
			}
			printf("displayMode:%d\n", displayMode);
			windowClearKey(0x44);
		}
	}
	return -1;
}
