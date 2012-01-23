#include <stdio.h>
#include <math.h>
#include <memory.h>

#include "ntscDecode.h"

typedef unsigned char Sample;  /* Sync = 4, Blank = 60, Black = 70, White = 200 */

#define SAMPLEDATA_MAX_SIZE (NTSC_SAMPLES_PER_LINE*NTSC_LINES_PER_FRAME)
#define INTERNAL_TEXTURE_MAX_SIZE (NTSC_SAMPLES_PER_LINE*NTSC_LINES_PER_FRAME)

static Sample s_sampleData[SAMPLEDATA_MAX_SIZE];

static unsigned int s_internalTexture[INTERNAL_TEXTURE_MAX_SIZE];
static unsigned int* s_outputTexture = NULL;
static int s_sampleAddIndex = 0;
static int s_sampleReadIndex = 0;

static Sample readerGetItem(const int offset)
{
	int ind = offset + s_sampleReadIndex;
	while (ind >= SAMPLEDATA_MAX_SIZE)
	{
		ind -= SAMPLEDATA_MAX_SIZE;
	}
	return s_sampleData[ind];
}

static void consumerRead(const int amount)
{
	s_sampleReadIndex += amount;
	while (s_sampleReadIndex >= SAMPLEDATA_MAX_SIZE)
	{
		s_sampleReadIndex -= SAMPLEDATA_MAX_SIZE;
	}
}

typedef struct DecodeNTSC
{
    unsigned char* _dataData;
    int _dataPitch;

    int _minSamplesPerLine;
    int _maxSamplesPerLine;
    int _minLinesPerField;
    int _maxLinesPerField;

    float _brightness;
    float _contrast;
    float _saturation;
    float _tint;
    float _horizontalSize;
    float _horizontalPosition;
    float _verticalSize;
    float _verticalPosition;
    int _verticalHold;
    int _horizontalHold;
    float _bloomFactor;
    float _scanLineHeight;

    float _linesVisible;
    float _lineTop;
    int _driftLines;
    int _driftSamples;
    float _active;
    float _preActive;
    int _colorBurstStart;

    float _baseLoad;
    float _crtLoad;

    float _colorBurstPhase[4];
    float _lockedColorBurstPhase[4];

    int _lastTime;
    int _frames;

    int _phase;

    int _line;
    int _foundVerticalSync;
    int _verticalSync;
    float _verticalSyncPhase;

    float _lefts[14+263];
    float _widths[14+263];

    int _delay[19+NTSC_SAMPLES_PER_LINE+8];
    int _gamma[256];
    int* _gamma0;

    int _topLine;
    int _bottomLine;

    int _linePeriod;

    int _hysteresisCount;
    int _colorMode;

		int _n;
} DecodeNTSC;

void decodeNTSCpostField(DecodeNTSC* const pDecodeNTSC)
{
	int lines = pDecodeNTSC->_bottomLine - pDecodeNTSC->_topLine;
	int y;
	for (y = 0; y < lines; ++y) 
	{
		int line = y + pDecodeNTSC->_topLine;
		float top = ((float)(line) - pDecodeNTSC->_verticalSyncPhase - pDecodeNTSC->_lineTop);
		int topI = (int)top;
		memcpy(s_outputTexture+(y+topI)*NTSC_SAMPLES_PER_LINE, s_internalTexture+y*NTSC_SAMPLES_PER_LINE*2, NTSC_SAMPLES_PER_LINE*sizeof(int));
		if (0)
		{
			/*
				 float left = _lefts[y]/static_cast<float>(_maxSamplesPerLine);
				 float right = left + _widths[y]/static_cast<float>(_maxSamplesPerLine);
				 */
		}
		memset(s_internalTexture+y*NTSC_SAMPLES_PER_LINE*2, 0, NTSC_SAMPLES_PER_LINE*sizeof(int));
	}
	pDecodeNTSC->_line = 0;
	pDecodeNTSC->_foundVerticalSync = 0;
	pDecodeNTSC->_crtLoad = pDecodeNTSC->_baseLoad;
	pDecodeNTSC->_verticalSync = 0;
}

/* Process a single scanline here */
static void decodeNTSCprocess(DecodeNTSC* const pDecodeNTSC, const int displayMode)
{
	/* Find the horizontal sync position. */
	int offset = 0;
	int i;
	int linePeriod;
	int p0;
	float total;
	float colorBurstGain;
	int phaseCorrelation;
	float colorBurstI0;
	float colorBurstQ0;
	float hf;
	int colorMode;
	for (i = 0; i < pDecodeNTSC->_driftSamples*2; ++i, ++offset)
	{
		if ((int)(readerGetItem(offset)) + (int)(readerGetItem(offset + 1)) < (int)pDecodeNTSC->_horizontalHold*2)
		{
			break;
		}
	}

	/* We use a phase-locked loop like real hardware does, in order to */
	/* avoid losing horizontal sync if the pulse is missing for a line or */
	/* two, and so that we get the correct "wobble" behavior. */
	linePeriod = pDecodeNTSC->_maxSamplesPerLine - offset;
	pDecodeNTSC->_linePeriod = (2*pDecodeNTSC->_linePeriod + linePeriod)/3;
	pDecodeNTSC->_linePeriod = clampInt(pDecodeNTSC->_minSamplesPerLine, pDecodeNTSC->_linePeriod, pDecodeNTSC->_maxSamplesPerLine);
	offset = pDecodeNTSC->_maxSamplesPerLine - pDecodeNTSC->_linePeriod;

	/* Find the vertical sync position. */
	if (!pDecodeNTSC->_foundVerticalSync)
	{
		int j;
		for (j = 0; j < pDecodeNTSC->_maxSamplesPerLine; j += 57) 
		{
			pDecodeNTSC->_verticalSync = ((pDecodeNTSC->_verticalSync*232)>>8) + (int)(readerGetItem(j)) - 60;
			if (pDecodeNTSC->_verticalSync < -pDecodeNTSC->_verticalHold || pDecodeNTSC->_line == 2*pDecodeNTSC->_driftLines) 
			{
				/* To render interlaced signals correctly, we need to */
				/* figure out where the vertical sync pulse happens */
				/* relative to the horizontal pulse. This determines the */
				/* vertical position of the raster relative to the screen. */
				pDecodeNTSC->_verticalSyncPhase = (float)(j)/(float)(pDecodeNTSC->_maxSamplesPerLine);
				/* Now we can find out which scanlines are at the top and */
				/* bottom of the screen. */
				pDecodeNTSC->_topLine = (int)(0.5f + pDecodeNTSC->_lineTop + pDecodeNTSC->_verticalSyncPhase);
				pDecodeNTSC->_bottomLine = (int)(1.5f + pDecodeNTSC->_linesVisible + pDecodeNTSC->_lineTop + pDecodeNTSC->_verticalSyncPhase);
				pDecodeNTSC->_line = 0;
				pDecodeNTSC->_foundVerticalSync = 1;
				break;
			}
		}
	}

	/* Determine the phase and strength of the color signal from the color */
	/* burst, which starts shortly after the horizontal sync pulse ends. */
	/* The color burst is 9 cycles long, and we look at the middle 5 cycles. */
	p0 = offset&~3;
	for (i = pDecodeNTSC->_colorBurstStart + 8; i < pDecodeNTSC->_colorBurstStart + 28; ++i) 
	{
		static const float colorBurstFadeConstant = 1.0f/128.0f;
		pDecodeNTSC->_colorBurstPhase[(i + pDecodeNTSC->_phase)&3] =
			pDecodeNTSC->_colorBurstPhase[(i + pDecodeNTSC->_phase)&3]*(1.0f - colorBurstFadeConstant) +
			(float)((int)(readerGetItem(p0 + i)) - 60)*colorBurstFadeConstant;
	}
	total = 0.1f;
	for (i = 0; i < 4; ++i)
	{
		total += pDecodeNTSC->_colorBurstPhase[i]*pDecodeNTSC->_colorBurstPhase[i];
	}

	colorBurstGain = 32.0f/sqrtf(total);
	phaseCorrelation = (offset + pDecodeNTSC->_phase)&3;
	colorBurstI0 = colorBurstGain*(pDecodeNTSC->_colorBurstPhase[2] - pDecodeNTSC->_colorBurstPhase[0])/16.0f;
	colorBurstQ0 = colorBurstGain*(pDecodeNTSC->_colorBurstPhase[3] - pDecodeNTSC->_colorBurstPhase[1])/16.0f;
	hf = colorBurstGain*(pDecodeNTSC->_colorBurstPhase[0] - pDecodeNTSC->_colorBurstPhase[1] + 
						pDecodeNTSC->_colorBurstPhase[2] - pDecodeNTSC->_colorBurstPhase[3]);
	colorMode = (colorBurstI0*colorBurstI0 + colorBurstQ0*colorBurstQ0) > 2.8 && hf < 16.0f;
	if (colorMode)
	{
		for (i = 0; i < 4; ++i)
		{
			pDecodeNTSC->_lockedColorBurstPhase[i] = pDecodeNTSC->_colorBurstPhase[i];
		}
		/* Color killer hysteresis: We only switch between colour mode and */
		/* monochrome mode if we stay in the new mode for 128 consecutive lines. */
		if (pDecodeNTSC->_colorMode != colorMode) 
		{
			pDecodeNTSC->_hysteresisCount++;
			if (pDecodeNTSC->_hysteresisCount == 128) 
			{
				pDecodeNTSC->_colorMode = colorMode;
				pDecodeNTSC->_hysteresisCount = 0;
			}
		}
	}
	else
	{
		pDecodeNTSC->_hysteresisCount = 0;
	}

	if (pDecodeNTSC->_foundVerticalSync && pDecodeNTSC->_line >= pDecodeNTSC->_topLine && pDecodeNTSC->_line < pDecodeNTSC->_bottomLine) 
	{
		float bloom;
		float horizontalSize;
		float samplesVisible;
		float sampleLeft;
		int start;
		int end;
		int brightness;
		unsigned int* destination;
		int yline;
		int totalSignal;

		yline = pDecodeNTSC->_line - pDecodeNTSC->_topLine;

		/* Lines with high amounts of brightness cause more load on the */
		/* horizontal oscillator which decreases horizontal deflection, */
		/* causing "blooming" (increase in width). */
		totalSignal = 0;
		for (i = 0; i < pDecodeNTSC->_active; ++i)
		{
			totalSignal += (int)(readerGetItem(offset + i)) - 60;
		}
		pDecodeNTSC->_crtLoad = 0.4f*pDecodeNTSC->_crtLoad + 0.6f*(pDecodeNTSC->_baseLoad + ((float)totalSignal - 42000.0f)/140000.0f);
		bloom = clampFloat(-2.0f, pDecodeNTSC->_bloomFactor*pDecodeNTSC->_crtLoad, 10.0f);
		horizontalSize = (1.0f - 6.3f*bloom/pDecodeNTSC->_active)*pDecodeNTSC->_horizontalSize;
		samplesVisible = pDecodeNTSC->_active*horizontalSize;
		sampleLeft = pDecodeNTSC->_preActive + pDecodeNTSC->_active*(0.5f + pDecodeNTSC->_horizontalPosition - horizontalSize/2.0f);
		pDecodeNTSC->_lefts[yline] = sampleLeft;
		pDecodeNTSC->_widths[yline] = samplesVisible;

		start = maxInt((int)(sampleLeft) - 10, 0);
		end = minInt((int)(sampleLeft + samplesVisible) + 10, pDecodeNTSC->_maxSamplesPerLine - offset);
		brightness = (int)(pDecodeNTSC->_brightness*100.0 - 7.5f*256.0f*pDecodeNTSC->_contrast)<<8;
		destination = (unsigned int*)(pDecodeNTSC->_dataData + yline*pDecodeNTSC->_dataPitch) + start;

		if (pDecodeNTSC->_colorMode) 
		{
			int x;
			int* p;
			int sp;
			int yContrast = (int)(pDecodeNTSC->_contrast*1463.0f);
			float radians = (float)(M_PI)/180.0f;
			float tintI = -cosf((103.0f + pDecodeNTSC->_tint)*radians);
			float tintQ = sinf((103.0f + pDecodeNTSC->_tint)*radians);
			int iqMultipliers[4];
			float colorBurstI = pDecodeNTSC->_lockedColorBurstPhase[(2 + phaseCorrelation)&3] - pDecodeNTSC->_lockedColorBurstPhase[(0 + phaseCorrelation)&3];
			float colorBurstQ = pDecodeNTSC->_lockedColorBurstPhase[(3 + phaseCorrelation)&3] - pDecodeNTSC->_lockedColorBurstPhase[(1 + phaseCorrelation)&3];

			iqMultipliers[0] = (int)((colorBurstI*tintI - colorBurstQ*tintQ)*pDecodeNTSC->_saturation*pDecodeNTSC->_contrast*colorBurstGain*0.352f);
			iqMultipliers[1] = (int)((colorBurstQ*tintI + colorBurstI*tintQ)*pDecodeNTSC->_saturation*pDecodeNTSC->_contrast*colorBurstGain*0.352f);
			iqMultipliers[2] = -iqMultipliers[0];
			iqMultipliers[3] = -iqMultipliers[1];
			p = &pDecodeNTSC->_delay[pDecodeNTSC->_maxSamplesPerLine];
			for (x = 0; x < 19; ++x)
			{
				p[x] = 0;
			}
			sp = offset + start;
			for (x = start; x < end; ++x, --p) 
			{
				int s;
				int C;
				int Y;
				int I;
				int Q;
				int R;
				int G;
				int B;
				unsigned int red = 0;
				unsigned int green = 0;
				unsigned int blue = 0;
				/* We use a low-pass Finite Impulse Response filter to */
				/* remove high frequencies (including the color carrier */
				/* frequency) from the signal. We could just keep a */
				/* 4-sample running average but that leads to sharp edges in the resulting image. */
				s = (int)(readerGetItem(sp++)) - 60;
				p[0] = s;
				Y = (p[6] + p[0] + ((p[5] + p[1])<<2) + 7*(p[4] + p[2]) + (p[3]<<3));
				Y = Y*yContrast + brightness;
				C = (s<<16) - Y;

				p[6] = s*iqMultipliers[x&3];
				I = p[12] + p[6] + ((p[11] + p[7])<<2) + 7*(p[10] + p[8]) + (p[9]<<3);

				p[12] = s*iqMultipliers[(x + 3)&3];
				Q = p[18] + p[12] + ((p[17] + p[13])<<2) + 7*(p[16] + p[14]) + (p[15]<<3);

				R = pDecodeNTSC->_gamma0[clampInt(0, (Y + 243*I + 160*Q)>>16, 255)];
				G = pDecodeNTSC->_gamma0[clampInt(0, (Y -  71*I - 164*Q)>>16, 255)];
				B = pDecodeNTSC->_gamma0[clampInt(0, (Y - 283*I + 443*Q)>>16, 255)];

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
					Y = clampInt(0, Y, 255<<16) >> 16;
					red = (unsigned int)Y;
					green = (unsigned int)Y;
					blue = (unsigned int)Y;
				}
				else if (displayMode == DISPLAY_CHROMA)
				{
					C = 128+(clampInt(-128<<16, C, 128<<16) >> 16);
					red = (unsigned int)C;
					green = (unsigned int)C;
					blue = (unsigned int)C;
				}
				else if (displayMode == DISPLAY_I)
				{
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
					int sampleValue = s + 60;
					red = (unsigned int)sampleValue;
					green = (unsigned int)sampleValue;
					blue = (unsigned int)sampleValue;
				}
				red = (unsigned int)clampInt(0, (int)red, 255);
				green = (unsigned int)clampInt(0, (int)green, 255);
				blue = (unsigned int)clampInt(0, (int)blue, 255);
				*(destination++) = (unsigned int)(0xff000000 | (red<<16) | (green<<8) | blue);
			}
		}
		else 
		{
			int x;
			int sp = offset + start;
			int yContrast = (int)(pDecodeNTSC->_contrast*46816.0f);
			unsigned int red = 0;
			unsigned int green = 0;
			unsigned int blue = 0;
			for (x = start; x < end; ++x) 
			{
				int s = (int)(readerGetItem(sp++)) - 60;
				int y = pDecodeNTSC->_gamma0[clampInt(0, (s*yContrast + brightness)>>16, 255)];
				if (displayMode == DISPLAY_RGB)
				{
					red = (unsigned int)y;
					green = (unsigned int)y;
					blue = (unsigned int)y;
				}
				else if (displayMode == DISPLAY_Y)
				{
					red = (unsigned int)y;
					green = (unsigned int)y;
					blue = (unsigned int)y;
				}
				else if (displayMode == DISPLAY_SIGNAL)
				{
					int sampleValue = s + 60;
					red = (unsigned int)sampleValue;
					green = (unsigned int)sampleValue;
					blue = (unsigned int)sampleValue;
				}
				red = (unsigned int)clampInt(0, (int)red, 255);
				green = (unsigned int)clampInt(0, (int)green, 255);
				blue = (unsigned int)clampInt(0, (int)blue, 255);
				*(destination++) = (unsigned int)(0xff000000 | (red<<16) | (green<<8) | blue);
			}
		}
	}
	offset += pDecodeNTSC->_minSamplesPerLine;
	pDecodeNTSC->_phase = (pDecodeNTSC->_phase + offset)&3;
	consumerRead(offset);

	++pDecodeNTSC->_line;
	if (pDecodeNTSC->_foundVerticalSync && pDecodeNTSC->_line == pDecodeNTSC->_minLinesPerField)
	{
		decodeNTSCpostField(pDecodeNTSC);
	}
}

static void decodeNTSCinit(DecodeNTSC* const pDecodeNTSC)
{
	int i;

	float samplesPerSecond;
	float us;
	float sync;
	float breezeway;
	float colorBurst;
	float backPorch;
	float frontPorch;
	float blanking;
	float line;

	float preSyncLines;
	float syncLines;
	float postSyncLines;
	float lines;
	float blankingLines;
	float activeLines;

	pDecodeNTSC->_phase=0;
	pDecodeNTSC->_foundVerticalSync=0;
	pDecodeNTSC->_line=0;
	pDecodeNTSC->_baseLoad=0.5;
	pDecodeNTSC->_verticalSync=0;
	pDecodeNTSC->_hysteresisCount=0;
	pDecodeNTSC->_colorMode=0;

	/* TODO: make these user-settable */
	pDecodeNTSC->_brightness=0.06f;
	pDecodeNTSC->_contrast=3.0f;
	pDecodeNTSC->_saturation=0.7f;
	pDecodeNTSC->_tint=18.0f;
	pDecodeNTSC->_horizontalSize=0.95f;
	pDecodeNTSC->_horizontalPosition=0;
	pDecodeNTSC->_verticalSize=0.93f;
	pDecodeNTSC->_verticalPosition=-0.01f;
	pDecodeNTSC->_verticalHold=280;
	pDecodeNTSC->_horizontalHold=25;
	pDecodeNTSC->_bloomFactor=10.0f;

	samplesPerSecond = 157500000.0f/11.0f;
	us = samplesPerSecond/1000000.0f;  /* samples per microsecond */

	/* Horizontal times in samples. */
	sync = 4.7f*us;
	breezeway = 0.6f*us;
	pDecodeNTSC->_colorBurstStart = (int)(sync + breezeway);
	colorBurst = 2.5f*us;
	backPorch = 1.6f*us;
	frontPorch = 1.5f*us;
	blanking = sync + breezeway + colorBurst + backPorch + frontPorch;
	line = 910.0f;
	pDecodeNTSC->_active = line - blanking;
	pDecodeNTSC->_preActive = blanking - frontPorch;
	/* The following parameter specifies how many samples early or late the */
	/* horizontal sync pulse can be and still be recognized (assuming good */
	/* signal fidelity). This sets the angle of the diagonal lines that */
	/* vertical lines become when horizontal sync is lost. */
	pDecodeNTSC->_driftSamples = 8;
	pDecodeNTSC->_minSamplesPerLine = (int)(line - (float)pDecodeNTSC->_driftSamples);
	pDecodeNTSC->_maxSamplesPerLine = (int)(line + (float)pDecodeNTSC->_driftSamples);
	/* We always consume a scanline at a time, and we won't be called to */
	/* process until we're sure we have a scanline. */
	pDecodeNTSC->_n = pDecodeNTSC->_maxSamplesPerLine;
	pDecodeNTSC->_linePeriod = (int)(line);

	/* Vertical times in lines. */
	preSyncLines = 3.0f;
	syncLines = 3.0f;
	postSyncLines = 14.0f;
	lines = 262.5f;
	blankingLines = preSyncLines + syncLines + postSyncLines;
	activeLines = lines - blankingLines;
	pDecodeNTSC->_linesVisible = activeLines*pDecodeNTSC->_verticalSize;
	pDecodeNTSC->_lineTop = postSyncLines + activeLines*(0.5f + pDecodeNTSC->_verticalPosition - pDecodeNTSC->_verticalSize/2.0f);
	/* The following parameter specifies how many lines early or late the */
	/* vertical sync pulse can be and still be recognized (assuming good */
	/* signal fidelity). This sets the "roll speed" of the picture when */
	/* vertical sync is lost. Empirically determined from video of an IBM 5153 monitor. */
	pDecodeNTSC->_driftLines = 14;
	pDecodeNTSC->_minLinesPerField = (int)(lines - (float)pDecodeNTSC->_driftLines);
	pDecodeNTSC->_maxLinesPerField = (int)(lines + (float)pDecodeNTSC->_driftLines);

	for (i = 0; i < 4; ++i)
		pDecodeNTSC->_colorBurstPhase[i] = pDecodeNTSC->_lockedColorBurstPhase[i] = 0;

	pDecodeNTSC->_crtLoad = pDecodeNTSC->_baseLoad;

	for (i = 0; i < 256; ++i)
		pDecodeNTSC->_gamma[i] = (int)(pow((float)(i)/255.0f, 1.9f)*255.0f);
	pDecodeNTSC->_gamma0 = &pDecodeNTSC->_gamma[0];
}


static DecodeNTSC s_decodeNTSC;

void crtSimInit(unsigned int* pVideoMemoryBGRA)
{
	decodeNTSCinit(&s_decodeNTSC);
	s_decodeNTSC._dataPitch = sizeof(int) * NTSC_SAMPLES_PER_LINE * 2;
	s_decodeNTSC._dataData = (unsigned char*)s_internalTexture;
	s_outputTexture = pVideoMemoryBGRA;
	s_sampleAddIndex = 0;
	s_sampleReadIndex = 0;
	memset(s_internalTexture, 0, sizeof(int)*INTERNAL_TEXTURE_MAX_SIZE);
}

void crtSimAddSample(const unsigned char sample)
{
	s_sampleData[s_sampleAddIndex] = sample;
	s_sampleAddIndex++;
	if (s_sampleAddIndex >= SAMPLEDATA_MAX_SIZE)
	{
		s_sampleAddIndex = 0;
	}
}

void crtSimTick(const int displayMode)
{
	do
	{	
  	decodeNTSCprocess(&s_decodeNTSC, displayMode);
	} while (s_decodeNTSC._line != 0 || s_decodeNTSC._foundVerticalSync);
}

