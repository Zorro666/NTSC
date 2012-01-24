#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ntscDecode.h"
#include "ntscDecodeCrtsim.h"

#define DECODE_JAKE			(0)
#define DECODE_CRTSIM		(1)

#define NTSC_COLOUR_CARRIER (3579545.0f)
#define NTSC_SAMPLE_RATE (NTSC_COLOUR_CARRIER*4.0f)

#define	NTSC_Y_LPF_CUTOFF (3.0f*1000.0f*1000.0f)

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

const char* const DISPLAY_FLAGS[] = {
	"INTERLACED"
	"INVALID"
	};

static float s_LPFY_inSignal[3];
static float s_LPFY_outY[3];
static float s_LPFY_a[3];
static float s_LPFY_b[2];
static float s_CHROMA_T = 0.0f;
static int s_displayModeFlags = DISPLAY_RGB | (DISPLAY_INTERLACED << 16);
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
	const float r = 1.000f;
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

static int s_yLPF[7];

static void decodeSignalY(const int compositeSignal, int* outY)
{
	int Y = 0;

	/* Y = LPF(signal) : 6MHz low-pass */
	s_LPFY_inSignal[0] = s_LPFY_inSignal[1];
	s_LPFY_inSignal[1] = s_LPFY_inSignal[2];
	s_LPFY_inSignal[2] = (float)compositeSignal;
	s_LPFY_outY[0] = s_LPFY_outY[1];
	s_LPFY_outY[1] = s_LPFY_outY[2];
	/* inSignal = n-2, n-1, n, outY = n-2, n-1, n */
	lowPass(s_LPFY_a, s_LPFY_b, s_LPFY_inSignal, s_LPFY_outY);
	Y = (int)s_LPFY_outY[2];

	s_yLPF[6] = s_yLPF[5];
	s_yLPF[5] = s_yLPF[4];
	s_yLPF[4] = s_yLPF[3];
	s_yLPF[3] = s_yLPF[2];
	s_yLPF[2] = s_yLPF[1];
	s_yLPF[1] = s_yLPF[0];
	s_yLPF[0] = compositeSignal;
	Y = (s_yLPF[6] + s_yLPF[0] + ((s_yLPF[5] + s_yLPF[1])<<2) + 7*(s_yLPF[4] + s_yLPF[2]) + (s_yLPF[3]<<3));
	Y = Y / 32;

	Y = clampInt(0, Y, 255);

	*outY = Y;
}

static int s_iLPF[7];
static int s_qLPF[7];
static int s_jakeI = 0;
static float s_jakeVals[4];

static void decodeSignalIQ(const int compositeSignal, int* outI, int* outQ)
{
	int I = 0;
	int Q = 0;
	float sinValue = 0;
	float cosValue = 0;
	float sinColourCarrier;
	float cosColourCarrier;
	float chromaValue = (float)compositeSignal;

	/* demodulate chroma to I, Q */
	const float COLOUR_CARRIER_DELTA_T = (float)(2.0f * M_PI * NTSC_COLOUR_CARRIER / NTSC_SAMPLE_RATE);
	sinColourCarrier = sinf(s_CHROMA_T);
	cosColourCarrier = -cosf(s_CHROMA_T);
	s_CHROMA_T += COLOUR_CARRIER_DELTA_T;

	/*
	printf("sin:%f cos:%f jakeI:%d %f %f\n", sinColourCarrier, cosColourCarrier, s_jakeI, s_jakeVals[s_jakeI&3], s_jakeVals[(s_jakeI+3)&3]);
	*/
	sinColourCarrier = s_jakeVals[s_jakeI&3];
	cosColourCarrier = s_jakeVals[(s_jakeI+3)&3];
	s_jakeI++;

	sinValue = chromaValue * sinColourCarrier;
	cosValue = chromaValue * cosColourCarrier;

	s_iLPF[6] = s_iLPF[5];
	s_iLPF[5] = s_iLPF[4];
	s_iLPF[4] = s_iLPF[3];
	s_iLPF[3] = s_iLPF[2];
	s_iLPF[2] = s_iLPF[1];
	s_iLPF[1] = s_iLPF[0];
	s_iLPF[0] = (int)sinValue * 32;
	I = (s_iLPF[6] + s_iLPF[0] + ((s_iLPF[5] + s_iLPF[1])<<2) + 7*(s_iLPF[4] + s_iLPF[2]) + (s_iLPF[3]<<3));

	s_qLPF[6] = s_qLPF[5];
	s_qLPF[5] = s_qLPF[4];
	s_qLPF[4] = s_qLPF[3];
	s_qLPF[3] = s_qLPF[2];
	s_qLPF[2] = s_qLPF[1];
	s_qLPF[1] = s_qLPF[0];
	s_qLPF[0] = (int)cosValue * 32;
	Q = (s_qLPF[6] + s_qLPF[0] + ((s_qLPF[5] + s_qLPF[1])<<2) + 7*(s_qLPF[4] + s_qLPF[2]) + (s_qLPF[3]<<3));

	*outI = I;
	*outQ = Q;
}

static int s_lineOdd = 0;
static void lineInit(void)
{
	int i;
	for (i = 0; i < 3; i++)
	{
		s_LPFY_inSignal[i] = 0.0f;
		s_LPFY_outY[i] = 0.0f;
	}
	for (i = 0; i < 7; i++)
	{
		s_yLPF[i] = 0;
		s_iLPF[i] = 0;
		s_qLPF[i] = 0;
	}
	s_CHROMA_T = 0.0f;
 	s_CHROMA_T += (float)(M_PI / 180.0f) * 33.0f;
	s_CHROMA_T = 0.0f;

	s_jakeI = 0;
	s_jakeVals[0] = -1.0f;
	s_jakeVals[1] = -0.0f;
	s_jakeVals[2] = +1.0f;
	s_jakeVals[3] = +0.0f;
	s_jakeI = 2*s_lineOdd;
	s_lineOdd ^= 1;
}


/* Public API */
void ntscDecodeInit(unsigned int* pVideoMemoryBGRA)
{
	int displayMode;
	int displayFlags;
	s_decodeOption = DECODE_CRTSIM;
	s_decodeOption = DECODE_JAKE;

	s_displayModeFlags = DISPLAY_RGB | (DISPLAY_INTERLACED << 16);
	s_pVideoMemoryBGRA = pVideoMemoryBGRA;
	s_pixelPos = 0;

	computeLowPassCoeffs(s_LPFY_a, s_LPFY_b, NTSC_Y_LPF_CUTOFF, NTSC_SAMPLE_RATE);
	printf("LPFY coeffs\n");
	printf("a[0]:%f a[1]:%f a[2]:%f b[0]:%f b[1]:%f\n", s_LPFY_a[0], s_LPFY_a[1], s_LPFY_a[2], s_LPFY_b[0], s_LPFY_b[1]);

	displayMode = s_displayModeFlags & 0xFFFF;
	displayFlags = s_displayModeFlags >> 16;
	printf("displayMode:'%s' (%d) displayFlags:0x%X\n", DISPLAY_MODES[displayMode], displayMode, displayFlags);

	crtSimInit(pVideoMemoryBGRA);
}

void ntscDecodeAddSample(const unsigned char sampleValue)
{
	if (s_decodeOption == DECODE_CRTSIM)
	{
		crtSimAddSample(sampleValue);
	}

	if (s_decodeOption == DECODE_JAKE)
	{
		const int displayMode = s_displayModeFlags & 0xFFFF;
		/*const int displayFlags = s_displayModeFlags >> 16;*/
		int pixelPos = s_pixelPos;
		int compositeSignal;

		int Y;
		int C;
		int I;
		int Q;
		int R;
		int G;
		int B;

		unsigned int red = 0;
		unsigned int green = 0;
		unsigned int blue = 0;
		unsigned int alpha = 0xFF;

		unsigned int* texture = s_pVideoMemoryBGRA;

		if (pixelPos == 0)
		{
			s_lineOdd = 0;
		}
		if (pixelPos%910 == 0)
		{
			lineInit();
		}

		compositeSignal = sampleValue - 60;
		decodeSignalY(compositeSignal, &Y);
		C = compositeSignal - Y;
		decodeSignalIQ(compositeSignal, &I, &Q);

		R = (int)((float)Y + 0.9563f * (float)I/256.0f + 0.6210f * (float)Q/256.0f);
		G = (int)((float)Y - 0.2721f * (float)I/256.0f - 0.6474f * (float)Q/256.0f);
		B = (int)((float)Y - 1.1070f * (float)I/256.0f + 1.7406f * (float)Q/256.0f);
		if (displayMode == DISPLAY_RGB)
		{
			red = (unsigned int)R;
			green = (unsigned int)G;
			blue = (unsigned int)B;
		}
		if (displayMode == DISPLAY_RED)
		{
			red = (unsigned int)R;
			green = 0;
			blue = 0;
		}
		else if (displayMode == DISPLAY_GREEN)
		{
			red = 0;
			green = (unsigned int)G;
			blue = 0;
		}
		else if (displayMode == DISPLAY_BLUE)
		{
			red = 0;
			green = 0;
			blue = (unsigned int)B;
		}
		else if (displayMode == DISPLAY_Y)
		{
			Y = clampInt(0, Y, 255<<0) >> 0;
			red = (unsigned int)Y;
			green = (unsigned int)Y;
			blue = (unsigned int)Y;
		}
		else if (displayMode == DISPLAY_CHROMA)
		{
			C = 128+(clampInt(-128<<1, C, 128<<1) >> 1);
			red = (unsigned int)C;
			green = (unsigned int)C;
			blue = (unsigned int)C;
		}
		else if (displayMode == DISPLAY_I)
		{
			/*printf("I:%d Q:%d\n", I, Q);*/
			I = 128+(clampInt(-128<<8, I, 128<<8) >> 8);
			red = (unsigned int)I;
			green = (unsigned int)I;
			blue = (unsigned int)I;
		}
		else if (displayMode == DISPLAY_Q)
		{
			Q = 128+(clampInt(-128<<8, Q, 128<<8) >> 8);
			red = (unsigned int)Q;
			green = (unsigned int)Q;
			blue = (unsigned int)Q;
		}
		else if (displayMode == DISPLAY_SIGNAL)
		{
			red = (unsigned int)sampleValue;
			green = (unsigned int)sampleValue;
			blue = (unsigned int)sampleValue;
		}
		red = (unsigned int)clampInt(0, (int)red, 255);
		green = (unsigned int)clampInt(0, (int)green, 255);
		blue = (unsigned int)clampInt(0, (int)blue, 255);
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
	int displayMode = s_displayModeFlags & 0xFFFF;
	int displayFlags = s_displayModeFlags >> 16;

	const int oldDisplayMode = displayMode;
	const int oldDisplayFlags = displayFlags;

	if (s_decodeOption == DECODE_CRTSIM)
	{
		crtSimTick(s_displayModeFlags);
		crtSimTick(s_displayModeFlags);
	}

	if (windowCheckKey('J'))
	{
		windowClearKey('J');
		if (s_decodeOption == DECODE_JAKE)
		{
			s_decodeOption = DECODE_CRTSIM;
		}
		else
		{
			s_decodeOption = DECODE_JAKE;
		}
	}
	if (windowCheckKey('0'))
	{
		windowClearKey('0');
		memset(s_pVideoMemoryBGRA, 0, 4 * NTSC_SAMPLES_PER_LINE * NTSC_LINES_PER_FRAME);
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
	if (windowCheckKey('D'))
	{
		windowClearKey('D');
		displayFlags ^= (DISPLAY_INTERLACED);
	}
	if (displayMode > DISPLAY_MAX)
	{
		displayMode = 0;
	}
	if (displayMode != oldDisplayMode)
	{
		printf("displayMode:'%s' (%d)\n", DISPLAY_MODES[displayMode], displayMode);
	}
	if (displayFlags != oldDisplayFlags)
	{
		printf("displayFlags:'%s' (%d)\n", DISPLAY_MODES[displayMode], displayFlags);
	}
	s_displayModeFlags = displayMode | (displayFlags << 16);
}

