#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ntscDecode.h"
#include "ntscDecodeCrtsim.h"

#define DECODE_JAKE			(0)
#define DECODE_CRTSIM		(1)

#define NTSC_COLOUR_CARRIER (3579545.0f)
#define NTSC_SAMPLE_RATE (NTSC_COLOUR_CARRIER*4.0f)

#define	NTSC_Y_LPF_CUTOFF (2.0f*1000.0f*1000.0f)

const char* const DISPLAY_MODES[] = {
	"RGB",
	"Y_SIGNAL",
	"CHROMA_SIGNAL",
	"I_SIGNAL",
	"Q_SIGNAL",
	"RAW SIGNAL",
	"RED CHANNEL",
	"GREEN CHANNEL",
	"BLUE CHANNEL",
	"INVALID"
	};

static float s_LPFY_inSignal[3];
static float s_LPFY_outY[3];
static float s_LPFY_a[3];
static float s_LPFY_b[2];
static float s_CHROMA_T = 0.0f;
static int s_displayMode = DISPLAY_RGB;
static unsigned int* s_pVideoMemoryBGRA = NULL;
static int s_pixelPos = 0;
static int s_decodeOption = DECODE_JAKE;

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

static void computeLowPassCoeffs(float a[3], float b[2], const float freq, const float sample_rate)
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
static void lowPass( const float a[3], const float b[2], const float* const pSamples, float* const pOutput)
{
	/* out(n) = a1 * in + a2 * in(n-1) + a3 * in(n-2) - b1*out(n-1) - b2*out(n-2) */
	pOutput[2] = a[0]*pSamples[2] + a[1]*pSamples[1] + a[2]*pSamples[0] - b[0]*pOutput[1] - b[1]*pOutput[0];
}

static void decodeCompositeSignalYC(const unsigned char compositeSignal, unsigned char* outY, unsigned char* outC)
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

	*outY = (unsigned char)(Y - 60);
	*outC = C;
}

static void decodeChromaSignalIQ(const unsigned char chromaSignal, unsigned char* outI, unsigned char* outQ)
{
	unsigned char I = 0;
	unsigned char Q = 0;
	float sinValue = 0;
	float cosValue = 0;
	float sinColourCarrier;
	float cosColourCarrier;
	float chromaValue = (float)chromaSignal - 60.0f;

	/* demodulate chroma to I, Q */
	const float COLOUR_CARRIER_DELTA_T = (float)(2.0f * M_PI * NTSC_COLOUR_CARRIER / NTSC_SAMPLE_RATE);
	sinColourCarrier = 2.0f * sinf(s_CHROMA_T+(float)(33.0f*M_PI/180.0f));
	cosColourCarrier = 2.0f * cosf(s_CHROMA_T+(float)(33.0f*M_PI/180.0f));
	s_CHROMA_T += COLOUR_CARRIER_DELTA_T;

	sinValue = chromaValue * sinColourCarrier;
	cosValue = chromaValue * cosColourCarrier;

	if (sinValue < 0.0f)
	{
		sinValue = 0.0f;
	}
	if (sinValue > 200.0f)
	{
		sinValue = 200.0f;
	}
	if (cosValue < 0.0f)
	{
		cosValue = 0.0f;
	}
	if (cosValue > 200.0f)
	{
		cosValue = 200.0f;
	}
	I = (unsigned char)sinValue;
	Q = (unsigned char)cosValue;

	*outI = I;
	*outQ = Q;
}

static void lineInit(void)
{
	int i;
	for (i = 0; i < 3; i++)
	{
		s_LPFY_inSignal[i] = 0.0f;
		s_LPFY_outY[i] = 0.0f;
	}
	s_CHROMA_T = 0.0f + (float)(2.0f * M_PI / 360.0f) * 33.0f;
}


/* Public API */
void ntscDecodeInit(unsigned int* pVideoMemoryBGRA)
{
	s_decodeOption = DECODE_JAKE;
	s_decodeOption = DECODE_CRTSIM;

	s_displayMode = DISPLAY_RGB;
	s_pVideoMemoryBGRA = pVideoMemoryBGRA;
	s_pixelPos = 0;

	computeLowPassCoeffs(s_LPFY_a, s_LPFY_b, NTSC_Y_LPF_CUTOFF, NTSC_SAMPLE_RATE);
	printf("LPFY coeffs\n");
	printf("a[0]:%f a[1]:%f a[2]:%f b[0]:%f b[1]:%f\n", s_LPFY_a[0], s_LPFY_a[1], s_LPFY_a[2], s_LPFY_b[0], s_LPFY_b[1]);

	printf("displayMode:'%s' (%d)\n", DISPLAY_MODES[s_displayMode], s_displayMode);

	crtSimInit(pVideoMemoryBGRA);
}

void ntscDecodeAddSample(const unsigned char sampleValue)
{
	int pixelPos = s_pixelPos;
	int displayMode = s_displayMode;
	unsigned char red = 0;
	unsigned char green = 0;
	unsigned char blue = 0;
	unsigned char alpha = 0xFF;

	unsigned char Y = 0;
	unsigned char C = 0;
	unsigned char I = 0;
	unsigned char Q = 0;

	unsigned int* texture = s_pVideoMemoryBGRA;

	if (s_decodeOption == DECODE_CRTSIM)
	{
		crtSimAddSample(sampleValue);
	}

	if (s_decodeOption == DECODE_JAKE)
	{
		if (pixelPos%910 == 0)
		{
			lineInit();
		}


		decodeCompositeSignalYC(sampleValue, &Y, &C);
		decodeChromaSignalIQ(C, &I, &Q);

		red = (unsigned char)((float)Y + 0.9563f * (float)I + 0.6210f * (float)Q);
		green = (unsigned char)((float)Y - 0.2721f * (float)I - 0.6474f * (float)Q);
		blue = (unsigned char)((float)Y - 1.1070f * (float)I + 1.7406f * (float)Q);
		if (0)
		{
			red = (unsigned char)(((int)red * 255) / 200);
			green = (unsigned char)(((int)green * 255) / 200);
			blue = (unsigned char)(((int)blue * 255) / 200);
		}
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
		else if (displayMode == DISPLAY_RED)
		{
			green = 0;
			blue = 0;
		}
		else if (displayMode == DISPLAY_GREEN)
		{
			red = 0;
			blue = 0;
		}
		else if (displayMode == DISPLAY_BLUE)
		{
			red = 0;
			green = 0;
		}
		/* BGRA format */
		texture[pixelPos] = (unsigned int)((alpha<<24) | (red<<16) | (green<<8) | blue);
		pixelPos++;
		if (pixelPos >= NTSC_SAMPLES_PER_LINE * NTSC_LINES_PER_FRAME)
		{
			pixelPos = 0;
		}
		s_pixelPos = pixelPos;
	}
}

void ntscDecodeTick(void)
{
	int displayMode = s_displayMode;
	const int oldDisplayMode = displayMode;

	if (s_decodeOption == DECODE_CRTSIM)
	{
		crtSimTick(s_displayMode);
		crtSimTick(s_displayMode);
	}

	if (windowCheckKey('C'))
	{
		windowClearKey('C');
		s_decodeOption = DECODE_CRTSIM;
	}
	if (windowCheckKey('J'))
	{
		windowClearKey('J');
		s_decodeOption = DECODE_JAKE;
	}
	if (windowCheckKey('0'))
	{
		windowClearKey('0');
		memset(s_pVideoMemoryBGRA, 0, 4 * NTSC_SAMPLES_PER_LINE * NTSC_LINES_PER_FRAME);
	}

	if (windowCheckKey('D'))
	{
		windowClearKey('D');
		displayMode++;
	}
	if (windowCheckKey('R'))
	{
		windowClearKey('R');
		displayMode = DISPLAY_RGB;
	}
	if (windowCheckKey('1'))
	{
		windowClearKey('1');
		displayMode = DISPLAY_RED;
	}
	if (windowCheckKey('2'))
	{
		windowClearKey('2');
		displayMode = DISPLAY_GREEN;
	}
	if (windowCheckKey('3'))
	{
		windowClearKey('3');
		displayMode = DISPLAY_BLUE;
	}
	if (windowCheckKey('S'))
	{
		windowClearKey('S');
		displayMode = DISPLAY_SIGNAL;
	}
	if (windowCheckKey('Y'))
	{
		windowClearKey('Y');
		displayMode = DISPLAY_Y;
	}
	if (windowCheckKey('C'))
	{
		windowClearKey('C');
		displayMode = DISPLAY_CHROMA;
	}
	if (windowCheckKey('I'))
	{
		windowClearKey('I');
		displayMode = DISPLAY_I;
	}
	if (windowCheckKey('Q'))
	{
		windowClearKey('Q');
		displayMode = DISPLAY_Q;
	}
	if (displayMode > DISPLAY_MAX)
	{
		displayMode = 0;
	}
	if (displayMode != oldDisplayMode)
	{
		printf("displayMode:'%s' (%d)\n", DISPLAY_MODES[displayMode], displayMode);
		s_displayMode = displayMode;
	}
}

