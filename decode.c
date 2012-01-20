#include <stdio.h>
#include <math.h>
#include <malloc.h>

#include "window.h"

#define DISPLAY_RGB 		(0)
#define DISPLAY_Y 			(1)
#define DISPLAY_CHROMA 	(2)
#define DISPLAY_I 			(3)
#define DISPLAY_Q 			(4)
#define DISPLAY_SIGNAL 	(5)
#define DISPLAY_TEST 		(6)
#define DISPLAY_MAX 		(6)

#define NTSC_COLOUR_CARRIER (3579545.0f)
#define NTSC_SAMPLE_RATE (NTSC_COLOUR_CARRIER*4.0f)
#define NTSC_SAMPLES_PER_LINE (910)
#define NTSC_LINES_PER_FIELD (262)
#define NTSC_FIELDS_PER_IMAGE (2)

#define	NTSC_Y_LPF_CUTOFF (3.0f*1000.0f*1000.0f)

const unsigned int MAIN_WINDOW = 0;

const unsigned int MAIN_WIDTH = NTSC_SAMPLES_PER_LINE;
const unsigned int MAIN_HEIGHT = NTSC_LINES_PER_FIELD*NTSC_FIELDS_PER_IMAGE;

const char* const DISPLAY_MODES[] = {
	"RGB",
	"Y_SIGNAL",
	"CHROMA_SIGNAL",
	"I_SIGNAL",
	"Q_SIGNAL",
	"RAW SIGNAL",
	"TEST PATTERN",
	"INVALID"
	};

/*
Y = LPF2(video)
Sin = gen sin ( color burst(input))
I = LPF3( (video − Y ) ∗ Sin )
Q = LPF3( (video − Y ) ∗ Sin(2 :) )
*/

#if 0
Notch filter - be even better to use this instead of LPF - notch at the colour sub-carrier frequency

Parameters:
0 =< freq =< samplerate/2
0 =< q < 1 (The higher, the narrower)

AlgoAlgo=double pi = 3.141592654;
double sqrt2 = sqrt(2.0);

double freq = 2050; /* Change! (zero & pole angle) */
double q = 0.4;     /* Change! (pole magnitude) */

double z1x = cos(2*pi*freq/samplerate);
double a0a2 = (1-q)*(1-q)/(2*(fabs(z1x)+1)) + q;
double a1 = -2*z1x*a0a2;
double b1 = -2*z1x*q;
double b2 = q*q;
double reg0, reg1, reg2;

unsigned int streamofs;
reg1 = 0;
reg2 = 0;

/* Main loop */
for (streamofs = 0; streamofs < streamsize; streamofs++)
{
  reg0 = a0a2 * ((double)fromstream[streamofs]
                 + fromstream[streamofs+2])
       + a1 * fromstream[streamofs+1]
       - b1 * reg1
       - b2 * reg2;

  reg2 = reg1;
  reg1 = reg0;

  int temp = reg0;

  /* Check clipping */
  if (temp > 32767) {
    temp = 32767;
  } else if (temp < -32768) temp = -32768;

  /* Store new value */
  tostream[streamofs] = temp;
#endif

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

void computeLowPassCoeffs(float a[3], float b[2], const float freq, const float sample_rate)
{
	const float r = 1.414f;
	const float c = 1.0f / tanf((float)(M_PI * freq/sample_rate));

	a[0] = 1.0f / (1.0f + r*c + c*c);
	a[1] = 2.0f * a[0];
	a[2] = a[0];
	b[0] = 2.0f * (1.0f - c*c) * a[0];
	b[1] = (1.0f - r*c + c*c) * a[0];
}

/* pSamples = input[n-2], pOutput = output[n-2] */
void lowPass( const float a[3], const float b[2], const float* const pSamples, float* const pOutput)
{
	/* out(n) = a1 * in + a2 * in(n-1) + a3 * in(n-2) - b1*out(n-1) - b2*out(n-2) */
	pOutput[2] = a[0]*pSamples[2] + a[1]*pSamples[1] + a[2]*pSamples[0] - b[0]*pOutput[1] - b[1]*pOutput[0];
}

int ntscLoadData(const char* const fileName, unsigned char** pNtscDataPtr, unsigned int* pNtscDataSize)
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

static float s_LPFY_inSignal[3];
static float s_LPFY_outY[3];
static float s_LPFY_a[3];
static float s_LPFY_b[2];
static float s_CHROMA_T = 0.0f;

void ntscDecodeInit(void)
{
	int i;
	for (i = 0; i < 3; i++)
	{
		s_LPFY_inSignal[i] = 0.0f;
		s_LPFY_outY[i] = 0.0f;
	}
	s_CHROMA_T = 0.0f;
}

void ntscDecodeCompositeSignalYC(const unsigned char compositeSignal, unsigned char* outY, unsigned char* outC)
{
	unsigned char Y = 0;
	unsigned char C = 0;

	/* Y = LPF(signal) : 6MHz low-pass */
	s_LPFY_inSignal[0] = s_LPFY_inSignal[1];
	s_LPFY_inSignal[1] = s_LPFY_inSignal[2];
	s_LPFY_inSignal[2] = compositeSignal;
	s_LPFY_outY[0] = s_LPFY_outY[1];
	s_LPFY_outY[1] = s_LPFY_outY[2];
	lowPass(s_LPFY_a, s_LPFY_b, s_LPFY_inSignal, s_LPFY_outY);
	if (s_LPFY_outY[2] < 0.0f)
	{
		Y = 0;
	}
	else if (s_LPFY_outY[2] > 200.0f)
	{
		Y = 200;
	}
	else
	{
		Y = (unsigned char)s_LPFY_outY[2];
	}

	/* subtract to get chroma */
	C = (unsigned char)(compositeSignal - Y);

	*outY = Y;
	*outC = C;
}

void ntscDecodeChromaSignalIQ(const unsigned char chromaSignal, unsigned char* outI, unsigned char* outQ)
{
	unsigned char I = 0;
	unsigned char Q = 0;
	unsigned int value = 0;
	float sinColourCarrier;
	float cosColourCarrier;

	/* demodulate chroma to I, Q */
	const float COLOUR_CARRIER_DELTA_T = (float)(2.0f * M_PI * NTSC_COLOUR_CARRIER / NTSC_SAMPLE_RATE);
	s_CHROMA_T += COLOUR_CARRIER_DELTA_T;
	sinColourCarrier = sinf(s_CHROMA_T);
	cosColourCarrier = cosf(s_CHROMA_T);

	value = (unsigned int)(chromaSignal * sinColourCarrier);
	I = (unsigned char)value;

	value = (unsigned int)(chromaSignal * cosColourCarrier);
	Q = (unsigned char)value;

	*outI = I;
	*outQ = Q;
}

#define NTSC_COLOUR_CARRIER (3579545.0f)
int main(int argc, char* argv[])
{
	int i;
	int displayMode = DISPLAY_Y;
	unsigned char* ntscDataPtr = NULL;
	unsigned int ntscDataSize = 0;
	unsigned int ntscFilenameIndex = 0;
	const char* ntscFilename = NULL;

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
	if (ntscLoadData(ntscFilename, &ntscDataPtr, &ntscDataSize))
	{
		fprintf(stderr, "ERROR loading NTSC data '%s'\n", ntscFilename);
		return -1;
	}

	windowSetup(MAIN_WINDOW, MAIN_WIDTH, MAIN_HEIGHT);

	computeLowPassCoeffs(s_LPFY_a, s_LPFY_b, NTSC_Y_LPF_CUTOFF, NTSC_SAMPLE_RATE);
	printf("LPFY coeffs\n");
	printf("a[0]:%f a[1]:%f a[2]:%f b[0]:%f b[1]:%f\n", s_LPFY_a[0], s_LPFY_a[1], s_LPFY_a[2], s_LPFY_b[0], s_LPFY_b[1]);

	printf("displayMode:'%s' (%d)\n", DISPLAY_MODES[displayMode], displayMode);
	while (1)
	{
		unsigned int x;
		unsigned int y;
		unsigned char* const texture = windowGetVideoMemoryBGRA(MAIN_WINDOW);
		unsigned int sample;

		ntscDecodeInit();
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

				unsigned char Y = 0;
				unsigned char C = 0;
				unsigned char I = 0;
				unsigned char Q = 0;
				ntscDecodeCompositeSignalYC(sampleValue, &Y, &C);
				ntscDecodeChromaSignalIQ(C, &I, &Q);

				if (displayMode == DISPLAY_SIGNAL)
				{
					red = sampleValue;
					green = sampleValue;
					blue = sampleValue;
				}
				else if (displayMode == DISPLAY_Y)
				{
					red = (unsigned char)Y;
					green = (unsigned char)Y;
					blue = (unsigned char)Y;
				}
				else if (displayMode == DISPLAY_CHROMA)
				{
					red = (unsigned char)C;
					green = (unsigned char)C;
					blue = (unsigned char)C;
				}
				else if (displayMode == DISPLAY_I)
				{
					red = (unsigned char)I;
					green = (unsigned char)I;
					blue = (unsigned char)I;
				}
				else if (displayMode == DISPLAY_Q)
				{
					red = (unsigned char)Q;
					green = (unsigned char)Q;
					blue = (unsigned char)Q;
				}
				else if (displayMode == DISPLAY_RGB)
				{
					red = (unsigned char)(Y + 0.9563f * I + 0.6210f * Q);
					green = (unsigned char)(Y - 0.2721f * I - 0.6474f * Q);
					blue = (unsigned char)(Y - 1.1070f * I + 1.7406f * Q);
				}
				else if (displayMode == DISPLAY_TEST)
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
		if (windowCheckKey('D'))
		{
			windowClearKey('D');
			displayMode++;
			if (displayMode > DISPLAY_MAX)
			{
				displayMode = 0;
			}
			printf("displayMode:'%s' (%d)\n", DISPLAY_MODES[displayMode], displayMode);
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
			if (ntscLoadData(ntscFilename, &ntscDataPtr, &ntscDataSize))
			{
				fprintf(stderr, "ERROR loading NTSC data '%s'\n", ntscFilename);
				return -1;
			}
		}
		if (windowCheckKey('R'))
		{
			windowClearKey('R');
			displayMode = DISPLAY_RGB;
			printf("displayMode:'%s' (%d)\n", DISPLAY_MODES[displayMode], displayMode);
		}
		if (windowCheckKey('S'))
		{
			windowClearKey('S');
			displayMode = DISPLAY_SIGNAL;
			printf("displayMode:'%s' (%d)\n", DISPLAY_MODES[displayMode], displayMode);
		}
		if (windowCheckKey('Y'))
		{
			windowClearKey('Y');
			displayMode = DISPLAY_Y;
			printf("displayMode:'%s' (%d)\n", DISPLAY_MODES[displayMode], displayMode);
		}
		if (windowCheckKey('C'))
		{
			windowClearKey('C');
			displayMode = DISPLAY_CHROMA;
			printf("displayMode:'%s' (%d)\n", DISPLAY_MODES[displayMode], displayMode);
		}
		if (windowCheckKey('I'))
		{
			windowClearKey('I');
			displayMode = DISPLAY_I;
			printf("displayMode:'%s' (%d)\n", DISPLAY_MODES[displayMode], displayMode);
		}
		if (windowCheckKey('Q'))
		{
			windowClearKey('Q');
			displayMode = DISPLAY_Q;
			printf("displayMode:'%s' (%d)\n", DISPLAY_MODES[displayMode], displayMode);
		}
	}
	return -1;
}
