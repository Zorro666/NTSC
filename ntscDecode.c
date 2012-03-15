#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ntscDecode.h"
#include "ntscDecodeCrtsim.h"

#define NTSC_GAMMA (2.2f)

/* Sample values (/4 compared to NTSC reference levels because 8-bit instead of 10-bit */
/* Sync = 4, Blank = 60, Black = 70, White = 200 */
#define NTSC_VALUE_SYNC		(4)
#define NTSC_VALUE_BLANK	(60)
#define NTSC_VALUE_BLACK	(70)
#define NTSC_VALUE_WHITE	(200)

#define NTSC_VALUE_COLOUR_BURST_MIN		(28)
#define NTSC_VALUE_COLOUR_BURST_MAX		(88)

#define NTSC_COLOUR_BURST_START_SAMPLE (12)
#define NTSC_COLOUR_BURST_LENGTH_SAMPLE (4*9)
#define NTSC_BLANKING_SAMPLES (NTSC_COLOUR_BURST_START_SAMPLE+NTSC_COLOUR_BURST_LENGTH_SAMPLE+8)

#define DECODE_JAKE			(0)
#define DECODE_CRTSIM		(1)

#define NTSC_COLOUR_CARRIER (3579545.0f)
#define NTSC_SAMPLE_RATE (NTSC_COLOUR_CARRIER*4.0f)

#define	NTSC_Y_LPF_CUTOFF (3.0f*1000.0f*1000.0f)

/*
	HSYNC
	0V = SYNC
The format of the horizontal sync pulse varies. In the 525-line NTSC system it is a 4.85 µs-long pulse at 0 V
*/

/*
	 VSYNC
The format of such a signal in 525-line NTSC is:
pre-equalizing pulses (6 to start scanning odd lines, 5 to start scanning even lines)
long-sync pulses (5 pulses)
post-equalizing pulses (5 to start scanning odd lines, 4 to start scanning even lines)
Each pre- or post- equalizing pulse consists in half a scan line of black signal: 2 µs at 0 V, followed by 30 µs at 0.3 V.
Each long sync pulse consists in an equalizing pulse with timings inverted: 30 µs at 0 V, followed by 2 µs at 0.3 V.
0.3V = BLANK
*/

/* DO WE NEED TO CONSIDER BACK & FRONT PORCH - for hsync detection */
/* BACK PORCH is where colour burst is */

/* COLOUR BURST */
/*
burst signal = sin(at+b)
	Trying to find b have N samples of sin(at+b) and sin(at)
	Each sample is 90 deg apart

sin(at+b) = sin(at)*cos(b)+cos(at)*sin(b)
cos(at+b) = cos(at)*cos(b)-sin(at)*sin(b)
A = sin(at+b)*sin(at) = sin(at)*sin(at)*cos(b)+sin(at)*cos(at)*sin(b)
B = cos(at+b)*cos(at) = cos(at)*cos(at)*cos(b)-sin(at)*cos(at)*sin(b)
A + B = cos(b)
A = colour_burst_sample[0] * sin_colour_carrier
B = colour_burst_sample[1] * cos_colour_carrier

C = sin(at+b)*cos(at) = cos(at)*sin(at)*cos(b)+cos(at)*cos(at)*sin(b)
D = cos(at+b)*sin(at) = sin(at)*cos(at)*cos(b)-sin(at)*sin(at)*sin(b)
C - D = sin(b)
C = colour_burst_sample[2] * cos_colour_carrier
D = colour_burst_sample[3] * sin_colour_carrier

*/

static int s_decodeOption = DECODE_JAKE;
static int s_displayModeFlags = DISPLAY_RGB | (DISPLAY_INTERLACED << 16);
static unsigned int* s_pVideoMemoryBGRA = NULL;

static float s_CHROMA_T = 0.0f;
static int s_pixelPos = 0;
static int s_xpos = 0;
static int s_ypos = 0;
static int s_colourBurstTotal = 0;
static float s_colourBurstSamples[4];
static float s_colourBurstAvg = 0.0f;

static int s_syncSamples = 0;
static int s_blankSamples = 0;

static int s_sampleCounter = 0;
static int s_samplesPerField = 0;
static int s_ntscSamples = 0;

static int s_fieldCounter = 0;

static int s_vsyncFound = 0;
static int s_hsyncFound = 0;
static int s_hsyncPosition = 0;
static int s_colourBurstFound = 0;

static float s_contrast = 1.00f;
static float s_brightness = 0.0f;

static float s_gammaValue = 0.5f;

static const char* const DISPLAY_MODES[] = {
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

static const char* const DISPLAY_FLAGS[] = {
	"INTERLACED",
	"LEE MODE",
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
*/

static float s_LPFY_inSignal[3];
static float s_LPFY_outY[3];
static float s_LPFY_a[3];
static float s_LPFY_b[2];

static float s_yLPF[7];
static float s_iLPF[7];
static float s_qLPF[7];

static int s_jakeI = 0;
static float s_jakeVals[4];

static int s_gamma[256];

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

#if 0
/* pSamples = input[n-2], pOutput = output[n-2] */
static void lowPass( const float a[3], const float b[2], const float* const pSamples, float* const pOutput)
{
	/* out(n) = a1 * in + a2 * in(n-1) + a3 * in(n-2) - b1*out(n-1) - b2*out(n-2) */
	pOutput[2] = a[0]*pSamples[2] + a[1]*pSamples[1] + a[2]*pSamples[0] - b[0]*pOutput[1] - b[1]*pOutput[0];
}
#endif

static void decodeSignalY(const int compositeSignal, float* outY)
{
	float Y = 0;

#if 0
	/* Y = LPF(signal) : 6MHz low-pass */
	s_LPFY_inSignal[0] = s_LPFY_inSignal[1];
	s_LPFY_inSignal[1] = s_LPFY_inSignal[2];
	s_LPFY_inSignal[2] = (float)compositeSignal;
	s_LPFY_outY[0] = s_LPFY_outY[1];
	s_LPFY_outY[1] = s_LPFY_outY[2];
	/* inSignal = n-2, n-1, n, outY = n-2, n-1, n */
	lowPass(s_LPFY_a, s_LPFY_b, s_LPFY_inSignal, s_LPFY_outY);
	Y = (int)s_LPFY_outY[2];
#endif

	s_yLPF[6] = s_yLPF[5];
	s_yLPF[5] = s_yLPF[4];
	s_yLPF[4] = s_yLPF[3];
	s_yLPF[3] = s_yLPF[2];
	s_yLPF[2] = s_yLPF[1];
	s_yLPF[1] = s_yLPF[0];
	s_yLPF[0] = (float)compositeSignal;
	Y = (float)(s_yLPF[6] + s_yLPF[0] + ((s_yLPF[5] + s_yLPF[1])*4.0f) + 7.0f*(s_yLPF[4] + s_yLPF[2]) + (s_yLPF[3]*8.0f));
	Y = Y / 32.0f;

	Y = clampFloat(Y, 0.0f, 200.0f);

	*outY = Y;
}

static void decodeSignalIQ(const int compositeSignal, float* outI, float* outQ)
{
	float I = 0;
	float Q = 0;
	float sinValue = 0;
	float cosValue = 0;
	float sinColourCarrier;
	float cosColourCarrier;
	const float IQscaling = 5.0f/100.0f;
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

	sinValue = chromaValue * sinColourCarrier;
	cosValue = chromaValue * cosColourCarrier;

	s_iLPF[6] = s_iLPF[5];
	s_iLPF[5] = s_iLPF[4];
	s_iLPF[4] = s_iLPF[3];
	s_iLPF[3] = s_iLPF[2];
	s_iLPF[2] = s_iLPF[1];
	s_iLPF[1] = s_iLPF[0];
	s_iLPF[0] = sinValue;
	I = (float)(s_iLPF[6] + s_iLPF[0] + ((s_iLPF[5] + s_iLPF[1])*4.0f) + 7.0f*(s_iLPF[4] + s_iLPF[2]) + (s_iLPF[3]*8.0f));

	s_qLPF[6] = s_qLPF[5];
	s_qLPF[5] = s_qLPF[4];
	s_qLPF[4] = s_qLPF[3];
	s_qLPF[3] = s_qLPF[2];
	s_qLPF[2] = s_qLPF[1];
	s_qLPF[1] = s_qLPF[0];
	s_qLPF[0] = cosValue;
	Q = (float)(s_qLPF[6] + s_qLPF[0] + ((s_qLPF[5] + s_qLPF[1])*4.0f) + 7.0f*(s_qLPF[4] + s_qLPF[2]) + (s_qLPF[3]*8.0f));

	I = I * IQscaling;
	Q = Q * IQscaling;

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
	for (i = 0; i < 7; i++)
	{
		s_yLPF[i] = 0.0f;
		s_iLPF[i] = 0.0f;
		s_qLPF[i] = 0.0f;
	}
	s_CHROMA_T = 0.0f;
 	s_CHROMA_T += (float)(M_PI / 180.0f) * 33.0f;
	s_CHROMA_T = 0.0f;
}

static void computeGammaTable(void)
{
	int i;
	float gammaVal = s_gammaValue;
	for (i = 0; i < 256; i++)
	{
		s_gamma[i] = (int)(powf((float)(i)/255.0f, gammaVal)*255.0f);
	}
}

/* SVD: A  = U x E x V* */
/* A = m x n */
/* U = m x m */
/* E = m x n */
/* V = n x n */
/* in: a = A */
/* out: u = U, q = E, v = V */
extern int svd(unsigned int m, unsigned int n, int withu, int withv, float eps, float tol, float **a, float *q, float **u, float **v);

static void matrixPrintf(float** mat, unsigned int numRows, unsigned int numCols, const char* const name)
{
	unsigned int i;
	for (i = 0; i < numRows; i++)
	{
		unsigned int j;
		for (j = 0; j < numCols; j++)
		{
			printf("%s[%d][%d] = %f ", name, i, j, mat[i][j]);
		}
		printf("\n");
	}
}

static float** matrixMalloc(const unsigned int numRows, const unsigned int numCols)
{
	unsigned int i;
	float** result = NULL;
	result = malloc(sizeof(float*)*numRows);
	for (i = 0; i < numRows; i++)
	{
		unsigned int j;
		result[i] = malloc(sizeof(float)*numCols);
		for (j = 0; j < numCols; j++)
		{
			result[i][j] = 0.0f;
		}
	}
	return result;
}

static void matrixMultiply(float** result, float** left, float** right, 
													 const unsigned int numRowsLeft, const unsigned int numColsLeft, 
													 const unsigned int numColsRight)
{
	unsigned int i;
	for (i = 0; i < numRowsLeft; i++)
	{
		unsigned int j;
		for (j = 0; j < numColsRight; j++)
		{
			unsigned int k;
			result[i][j] = 0.0f;
			for (k = 0; k < numColsLeft; k++)
			{
				result[i][j] += left[i][k] * right[k][j];
			}
		}
	}
}

void computeColourBurstMatrices(void)
{
	const unsigned int M = 4;
	const unsigned int N = 2;

	float** UxE;
	float** E;
	float** jakeTemp1;

	float temp;

	unsigned int i;
	unsigned int j;
	int result;

	float* w;
	float** v;
	float** a;
	float** u;

	a = matrixMalloc(M, M);
	u = matrixMalloc(M, M);
	v = matrixMalloc(N, N);
	w = malloc(sizeof(float)*N);

	UxE = matrixMalloc(M, M);
	E = matrixMalloc(M, M);
	jakeTemp1 = matrixMalloc(M, M);
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < M; j++)
		{
			UxE[i][j] = 0.0f;
			E[i][j] = 0.0f;
			jakeTemp1[i][j] = UxE[i][j];
		}
	}

	for (i = 0; i < M; i++)
	{
		const float angle = (33.0f + ((float)(i+0)*90.0f)) * (float)M_PI/180.0f;
		const float sinC = (float)sinf(angle);
		const float cosC = (float)cosf(angle);
		a[i][0] = sinC;
		a[i][1] = cosC;
	}
	/*
	a[0][0] = 2.0f;
	a[0][1] = 1.0f;
	a[0][2] = 0.0f;
	a[0][3] = 0.0f;
	a[1][0] = 4.0f;
	a[1][1] = 3.0f;
	a[1][2] = 0.0f;
	a[1][3] = 0.0f;
	a[2][0] = 0.0f;
	a[2][1] = 0.0f;
	a[2][2] = 0.0f;
	a[2][3] = 0.0f;
	a[3][0] = 0.0f;
	a[3][1] = 0.0f;
	a[3][2] = 0.0f;
	a[3][3] = 0.0f;
	*/

	printf("Input\n");
	matrixPrintf(a, M, N, "a");

	/* SVD: A  = U x E x V* */
	/* out: a = U, w = E (just the non-zero values), v = V */
	/* A = m x n */
	/* U = m x m */
	/* E = m x n : 0 except on diagonal */
	/* V = n x n */
	result = svd(M, N, 1, 1, 1.0e-10f, 1.0e-10f, a, w, u, v);

	printf("SVD result = %d\n", result);

	matrixPrintf(u, M, M, "u");
	for (j = 0; j < N; j++)
	{
		printf("w[%d] = %f\n", j, w[j]);
	}
	matrixPrintf(v, N, N, "v");

	/* E = m x n : 0 except on diagonal */
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			float Evalue = 0.0f;
			if (i == j)
			{
				Evalue = w[i];
			}
			E[i][j] = Evalue;
		}
	}
	matrixPrintf(E, M, N, "E");

	/* UxE = U x E : U = m x m, E = m x n */
	matrixMultiply(UxE, u, E, M, N, N);
	matrixPrintf(UxE, M, N, "UxE");
	/* jakeTemp1 = (U x E) x V* = UxE * v : UxE = m x n, V* = n x n */
	temp = v[0][1];
	v[0][1] = v[1][0];
	v[1][0] = temp;
	matrixMultiply(jakeTemp1, UxE, v, M, N, N);
	matrixPrintf(jakeTemp1, M, N, "UxExV*");
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
	s_xpos = 0;
	s_ypos = 0;
	s_pixelPos = 0;

	computeColourBurstMatrices();

	computeLowPassCoeffs(s_LPFY_a, s_LPFY_b, NTSC_Y_LPF_CUTOFF, NTSC_SAMPLE_RATE);
	printf("LPFY coeffs\n");
	printf("a[0]:%f a[1]:%f a[2]:%f b[0]:%f b[1]:%f\n", s_LPFY_a[0], s_LPFY_a[1], s_LPFY_a[2], s_LPFY_b[0], s_LPFY_b[1]);

	displayMode = s_displayModeFlags & 0xFFFF;
	displayFlags = s_displayModeFlags >> 16;
	printf("displayMode:'%s' (%d) displayFlags:0x%X\n", DISPLAY_MODES[displayMode], displayMode, displayFlags);

	crtSimInit(pVideoMemoryBGRA);

	s_jakeI = 0;
	s_jakeVals[0] = +0.0f;
	s_jakeVals[1] = +0.0f;
	s_jakeVals[2] = -0.0f;
	s_jakeVals[3] = -0.0f;

	lineInit();

	s_ntscSamples = 0;

	s_syncSamples = 0;
	s_blankSamples = 0;

	s_sampleCounter = 0;

	s_fieldCounter = 0;
	s_samplesPerField = 0;

	s_vsyncFound = 0;
	s_hsyncFound = 0;
	s_hsyncPosition = 0;
	s_colourBurstFound = 0;

	computeGammaTable();
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
		const int displayFlags = s_displayModeFlags >> 16;
		int pixelPos = s_pixelPos;
		int compositeSignal;
		int syncFound = 0;
		int blankFound = 0;
		int blankSignal =0;
		int hsyncFound = 0;
		int vsyncFound = 0;

		int Y = 0;
		int C = 0;
		int I = 0;
		int Q = 0;
		int R = 0;
		int G = 0;
		int B = 0;

		unsigned int red = 0;
		unsigned int green = 0;
		unsigned int blue = 0;
		unsigned int alpha = 0xFF;

		unsigned int* texture = s_pVideoMemoryBGRA;

		s_ntscSamples++;
		s_sampleCounter++;
		s_samplesPerField++;
		if (s_samplesPerField > 3000)
		{
			s_vsyncFound = 0;
		}
		if (s_sampleCounter > 200)
		{
			s_hsyncFound = 0;
		}

		compositeSignal = sampleValue - NTSC_VALUE_BLANK;
		/* HSYNC is 4.85us long which is 69.5 NTSC samples */
		if (sampleValue <= NTSC_VALUE_SYNC)
		{
			syncFound = 1;
			s_syncSamples++;
		}
		if (sampleValue <= NTSC_VALUE_BLANK)
		{
			blankFound = 1;
			s_blankSamples++;
		}

		if ((s_syncSamples > 60) && (syncFound==0))
		{
			hsyncFound = 1;
		}
		if ((s_syncSamples > 380) && (syncFound==0))
		{
			vsyncFound = 1;
		}
		/* Colour burst is 2.5us long which is 36 NTSC samples (9 waves) */
		/* It is about ~1.0us after the HSYNC pulse ends which is 11 samples */
		if ((s_hsyncFound == 1) && (s_colourBurstFound == 0))
		{
			static int sampleAtStart = 0;
			static int sampleAtStart2 = 0;
			const int colourBurstStart = NTSC_COLOUR_BURST_START_SAMPLE;
			const int colourBurstEnd = colourBurstStart + NTSC_COLOUR_BURST_LENGTH_SAMPLE;
			const int colourBurstLookStart = colourBurstStart + 8;
			const int colourBurstLookEnd = colourBurstEnd - 8;
			if (s_sampleCounter < colourBurstLookStart)
			{
				s_colourBurstTotal = 0;
				s_colourBurstAvg = 0.0f;
				s_colourBurstSamples[0] = 0.0f;
				s_colourBurstSamples[1] = 0.0f;
				s_colourBurstSamples[2] = 0.0f;
				s_colourBurstSamples[3] = 0.0f;
				sampleAtStart = s_ntscSamples;
				sampleAtStart2 = s_sampleCounter;
			}
			if ((s_sampleCounter >= colourBurstLookStart) && (s_sampleCounter < colourBurstLookEnd))
			{
				int burstSampleIndex = (s_samplesPerField ) & 0x3;
				burstSampleIndex = (s_hsyncPosition + s_sampleCounter) & 0x3;

				burstSampleIndex = (s_sampleCounter - colourBurstLookStart);
				burstSampleIndex += 0*s_hsyncPosition;
				burstSampleIndex += sampleAtStart;
				burstSampleIndex += 0*sampleAtStart2;

				burstSampleIndex &= 0x3;
				/* We are in the colour burst phase */
				s_colourBurstTotal += (compositeSignal*compositeSignal);
				s_colourBurstSamples[burstSampleIndex] = s_colourBurstSamples[burstSampleIndex] * 0.5f + 0.5f*(float)compositeSignal;
			}
			if (s_sampleCounter >= colourBurstLookEnd)
			{
				if (s_colourBurstTotal > 1)
				{
					float bI;
					float bQ;
					float sinC = (float)sinf(33.0f * (float)M_PI/180.0f);
					float cosC = (float)cosf(33.0f * (float)M_PI/180.0f);
					const int numSamples = (colourBurstLookEnd-colourBurstLookStart);
					s_colourBurstSamples[0] /= (float)(numSamples/4);
					s_colourBurstSamples[1] /= (float)(numSamples/4);
					s_colourBurstSamples[2] /= (float)(numSamples/4);
					s_colourBurstSamples[3] /= (float)(numSamples/4);
					s_colourBurstAvg = (float)sqrtf((float)s_colourBurstTotal) / (float)(numSamples);

					bI = s_colourBurstSamples[2]-s_colourBurstSamples[0]; 
					bQ = s_colourBurstSamples[3]-s_colourBurstSamples[1];
					bI /= s_colourBurstAvg;
					bQ /= s_colourBurstAvg;
					/*
				 	printf("y:%d colourBurstTotal:%d avg:%f\n", s_ypos, s_colourBurstTotal, s_colourBurstAvg);
					*/
/*
C = sin(at+b)*cos(at) = cos(at)*sin(at)*cos(b)+cos(at)*cos(at)*sin(b)
D = cos(at+b)*sin(at) = sin(at)*cos(at)*cos(b)-sin(at)*sin(at)*sin(b)
C - D = sin(b)
C = colour_burst_sample[2] * cos_colour_carrier
D = colour_burst_sample[3] * sin_colour_carrier
bI = sin(at+b)
bQ = cos(at+b)

the colour burst looks to be sin(-33 deg + omega*t)
colour burst is 33 deg relative to Q

A good glossary and information place:
http://techpubs.sgi.com/library/dynaweb_docs/0530/SGI_Developer/books/Ind2Vid_PG/sgi_html/go01.html

cos(at+b) = (sin(180+at+b)-sin(at+b))*sin_colour_carrier - (cos(at+b+180)-cos(at+b))*cos_colour_carrier
cos(at+b) = sin(180+at+b)*sin_colour_carrier - cos(at+b+180)*cos_colour_carrier -
          = sin(at+b)*sin_colour_carrier - cos(at+b)*cos_colour_carrier
*/
					/*
					printf("jI:%d h:%d sam:%d sam2:%d y:%d n:%d a:%f burstSamples:%.3ff, %.3ff, %.3ff, %.3f\n", 
							s_jakeI&0x3, s_hsyncPosition&0x3, sampleAtStart&0x3, sampleAtStart2&0x3,
							s_ypos, 
							numSamples,
							s_colourBurstAvg,
							s_colourBurstSamples[0], s_colourBurstSamples[1], s_colourBurstSamples[2], s_colourBurstSamples[3]
							);
					*/
					/*
					printf("sinC:%f cosC:%f bI:%f bQ:%f vals:%f, %f\n", 
							sinC, cosC, bI, bQ, 
							(bI * cosC + bQ * sinC),
							(bI * sinC - bQ * cosC));
					*/

					s_jakeVals[0] = bI * sinC - bQ * cosC;
					s_jakeVals[1] = bI * cosC + bQ * sinC;
					s_jakeVals[2] = -s_jakeVals[0];
					s_jakeVals[3] = -s_jakeVals[1];
				}
				else
				{
					printf("Monochrome line - no colour burst y:%d\n", s_ypos);
				}
				s_colourBurstFound = 1;
			}
		}
		if (s_sampleCounter < NTSC_BLANKING_SAMPLES)
		{
			blankSignal = 1;
		}

		if (hsyncFound == 1)
		{
			/* HSYNC */
			/*printf("\tHSYNC x:%d y:%d\n", s_xpos, s_ypos);*/

			s_xpos = 0;
			lineInit();

			s_hsyncPosition = s_samplesPerField;
			s_sampleCounter = 0;
			if (displayFlags & DISPLAY_INTERLACED)
			{
				s_ypos += 2;
			}
			else
			{
				s_ypos++;
			}
			if (s_ypos >= NTSC_LINES_PER_FRAME)
			{
				s_ypos = NTSC_LINES_PER_FRAME-1;
			}
			pixelPos = s_ypos * NTSC_SAMPLES_PER_LINE;
			s_hsyncFound = 1;
			s_colourBurstFound = 0;
		}
		if ((vsyncFound == 1) && (s_vsyncFound == 0))
		{
			/* VSYNC */
			/*
			printf("VSYNC x:%d y:%d syncSamples:%d samples:%d\n", s_xpos, s_ypos, s_syncSamples, s_samplesPerField);
			*/
			s_fieldCounter++;
			s_ypos = s_fieldCounter & 0x1;
			if (displayFlags & DISPLAY_INTERLACED)
			{
			}
			else
			{
				s_ypos *= 262;
			}
			pixelPos = s_ypos * NTSC_SAMPLES_PER_LINE;
			/*
			printf("VSYNC newy:%d\n", s_ypos);
			*/
			s_samplesPerField = 0;
			s_vsyncFound = 1;
			if (displayFlags & DISPLAY_LEE_MODE)
			{
				/* For NES - for Lee */
				s_jakeI = s_fieldCounter & 0x1 ? 1 : 1;
			}
			else
			{
				/* For NTSC saved files */
				/*s_jakeI = s_fieldCounter & 0x1 ? 0 : 1;*/
			}
		}

		if (s_blankSamples >= 4)
		{
			blankSignal = 1;
		}

		if (blankSignal == 0)
		{
			float yval;
			float ival = 0.0f;
			float qval = 0.0f;
			decodeSignalY(compositeSignal, &yval);
			yval = yval * (255.0f/200.0f);
			yval = yval * s_contrast;
			yval = yval + s_brightness;
			if (s_colourBurstTotal > 1)
			{
				C = compositeSignal - (int)yval;
				decodeSignalIQ(compositeSignal, &ival, &qval);
			}

			R = (int)((float)yval + 0.9563f * ival + 0.6210f * qval);
			G = (int)((float)yval - 0.2721f * ival - 0.6474f * qval);
			B = (int)((float)yval - 1.1070f * ival + 1.7406f * qval);
			R = s_gamma[clampInt(R, 0, 255)];
			G = s_gamma[clampInt(G, 0, 255)];
			B = s_gamma[clampInt(B, 0, 255)];

			Y = (int)yval;
			I = (int)ival;
			Q = (int)qval;
		}
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
			Y = clampInt(Y, 0, 255<<0) >> 0;
			red = (unsigned int)Y;
			green = (unsigned int)Y;
			blue = (unsigned int)Y;
		}
		else if (displayMode == DISPLAY_CHROMA)
		{
			C = 128+(clampInt(C, -128<<1, 128<<1) >> 1);
			red = (unsigned int)C;
			green = (unsigned int)C;
			blue = (unsigned int)C;
		}
		else if (displayMode == DISPLAY_I)
		{
			/*printf("I:%d Q:%d\n", I, Q);*/
			I = 128+(clampInt(I, -128<<1, 128<<1) >> 1);
			red = (unsigned int)I;
			green = (unsigned int)I;
			blue = (unsigned int)I;
		}
		else if (displayMode == DISPLAY_Q)
		{
			Q = 128+(clampInt(Q, -128<<1, 128<<1) >> 1);
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
		red = (unsigned int)clampInt((int)red, 0, 255);
		green = (unsigned int)clampInt((int)green, 0, 255);
		blue = (unsigned int)clampInt((int)blue, 0, 255);

		if (((s_ypos == 200) || (s_ypos == 201)) && (s_xpos == 150))
		{
			if ((s_ypos > 0) && (s_xpos > 0))
			{
				/*
				printf("x:%d y:%d RGB:%d, %d, %d Value:%d Y:%d C:%d I:%d Q:%d\n", 
						s_xpos, s_ypos, red, green, blue, sampleValue, Y, C, I, Q);
				*/
			}
		}
		/* BGRA format */
		texture[pixelPos] = (unsigned int)((alpha<<24) | (red<<16) | (green<<8) | blue);

		if (syncFound == 0)
		{
			s_syncSamples = 0;
		}
		if (blankFound == 0)
		{
			s_blankSamples = 0;
		}

		s_xpos++;
		s_jakeI++;
		pixelPos++;
		if (s_xpos >= NTSC_SAMPLES_PER_LINE)
		{
			s_xpos = NTSC_SAMPLES_PER_LINE-1;
			pixelPos--;
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
		/*crtSimTick(s_displayModeFlags);*/
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
	if (windowCheckKey('L'))
	{
		windowClearKey('L');
		displayFlags ^= (DISPLAY_LEE_MODE);
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

