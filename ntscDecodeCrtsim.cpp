#include <stdio.h>
#include <math.h>
#include <memory.h>

typedef unsigned char Sample;  /* Sync = 4, Blank = 60, Black = 70, White = 200 */

#define NTSC_SAMPLES_PER_LINE (910)
#define NTSC_LINES_PER_FIELD (525)
#define SAMPLEDATA_MAX_SIZE (910*525)
#define INTERNAL_TEXTURE_MAX_SIZE (NTSC_SAMPLES_PER_LINE*NTSC_LINES_PER_FIELD)

static Sample s_sampleData[SAMPLEDATA_MAX_SIZE];

static unsigned int s_internalTexture[INTERNAL_TEXTURE_MAX_SIZE];
static unsigned int* s_outputTexture = NULL;
static int s_sampleAddIndex = 0;
static int s_sampleReadIndex = 0;

static Sample readerGetItem(const unsigned int offset)
{
	int index = offset + s_sampleReadIndex;
	while (index >= SAMPLEDATA_MAX_SIZE)
	{
		index -= SAMPLEDATA_MAX_SIZE;
	}
	return s_sampleData[index];
}

static void consumerRead(const int amount)
{
	s_sampleReadIndex += amount;
	while (s_sampleReadIndex >= SAMPLEDATA_MAX_SIZE)
	{
		s_sampleReadIndex -= SAMPLEDATA_MAX_SIZE;
	}
}

static float min(float a, float b)
{
	if (a < b)
	{
		return a;
	}
	return b;
}

static float max(float a, float b)
{
	if (a > b)
	{
		return a;
	}
	return b;
}

static float clamp(float low, float value, float high) 
{ 
	return min(max(low, value), high); 
}

static int min(int a, int b)
{
	if (a < b)
	{
		return a;
	}
	return b;
}

static int max(int a, int b)
{
	if (a > b)
	{
		return a;
	}
	return b;
}

static int clamp(int low, int value, int high) 
{ 
	return min(max(low, value), high); 
}

struct DecodeNTSC
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
    bool _foundVerticalSync;
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
    bool _colorMode;

		int _n;
};

void decodeNTSCpostField(DecodeNTSC* const pDecodeNTSC)
{
	int lines = pDecodeNTSC->_bottomLine - pDecodeNTSC->_topLine;
	for (int y = 0; y < lines; ++y) 
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
	}
	pDecodeNTSC->_line = 0;
	pDecodeNTSC->_foundVerticalSync = false;
	pDecodeNTSC->_crtLoad = pDecodeNTSC->_baseLoad;
	pDecodeNTSC->_verticalSync = 0;
}

/* Process a single scanline here */
static void decodeNTSCprocess(DecodeNTSC* const pDecodeNTSC)
{
	// Find the horizontal sync position.
	int offset = 0;
	for (int i = 0; i < pDecodeNTSC->_driftSamples*2; ++i, ++offset)
		if ((int)(readerGetItem(offset)) + (int)(readerGetItem(offset + 1)) < pDecodeNTSC->_horizontalHold*2)
			break;

	// We use a phase-locked loop like real hardware does, in order to
	// avoid losing horizontal sync if the pulse is missing for a line or
	// two, and so that we get the correct "wobble" behavior.
	int linePeriod = pDecodeNTSC->_maxSamplesPerLine - offset;
	pDecodeNTSC->_linePeriod = (2*pDecodeNTSC->_linePeriod + linePeriod)/3;
	pDecodeNTSC->_linePeriod = clamp(pDecodeNTSC->_minSamplesPerLine, pDecodeNTSC->_linePeriod, pDecodeNTSC->_maxSamplesPerLine);
	offset = pDecodeNTSC->_maxSamplesPerLine - pDecodeNTSC->_linePeriod;

	// Find the vertical sync position.
	if (!pDecodeNTSC->_foundVerticalSync)
		for (int j = 0; j < pDecodeNTSC->_maxSamplesPerLine; j += 57) {
			pDecodeNTSC->_verticalSync = ((pDecodeNTSC->_verticalSync*232)>>8) + (int)(readerGetItem(j)) - 60;
			if (pDecodeNTSC->_verticalSync < -pDecodeNTSC->_verticalHold || pDecodeNTSC->_line == 2*pDecodeNTSC->_driftLines) {
				// To render interlaced signals correctly, we need to
				// figure out where the vertical sync pulse happens
				// relative to the horizontal pulse. This determines the
				// vertical position of the raster relative to the screen.
				pDecodeNTSC->_verticalSyncPhase = (float)(j)/(float)(pDecodeNTSC->_maxSamplesPerLine);
				// Now we can find out which scanlines are at the top and
				// bottom of the screen.
				pDecodeNTSC->_topLine = (int)(0.5f + pDecodeNTSC->_lineTop + pDecodeNTSC->_verticalSyncPhase);
				pDecodeNTSC->_bottomLine = (int)(1.5f + pDecodeNTSC->_linesVisible + pDecodeNTSC->_lineTop + pDecodeNTSC->_verticalSyncPhase);
				pDecodeNTSC->_line = 0;
				pDecodeNTSC->_foundVerticalSync = true;
				break;
			}
		}

	// Determine the phase and strength of the color signal from the color
	// burst, which starts shortly after the horizontal sync pulse ends.
	// The color burst is 9 cycles long, and we look at the middle 5
	// cycles.
	int p0 = offset&~3;
	for (int i = pDecodeNTSC->_colorBurstStart + 8; i < pDecodeNTSC->_colorBurstStart + 28; ++i) {
		static const float colorBurstFadeConstant = 1.0f/128.0f;
		pDecodeNTSC->_colorBurstPhase[(i + pDecodeNTSC->_phase)&3] =
			pDecodeNTSC->_colorBurstPhase[(i + pDecodeNTSC->_phase)&3]*(1.0f - colorBurstFadeConstant) +
			(float)((int)(readerGetItem(p0 + i)) - 60)*colorBurstFadeConstant;
	}
	float total = 0.1f;
	for (int i = 0; i < 4; ++i)
		total += pDecodeNTSC->_colorBurstPhase[i]*pDecodeNTSC->_colorBurstPhase[i];
	float colorBurstGain = 32.0f/sqrtf(total);
	int phaseCorrelation = (offset + pDecodeNTSC->_phase)&3;
	float colorBurstI0 = colorBurstGain*(pDecodeNTSC->_colorBurstPhase[2] - pDecodeNTSC->_colorBurstPhase[0])/16.0f;
	float colorBurstQ0 = colorBurstGain*(pDecodeNTSC->_colorBurstPhase[3] - pDecodeNTSC->_colorBurstPhase[1])/16.0f;
	float hf = colorBurstGain*(pDecodeNTSC->_colorBurstPhase[0] - pDecodeNTSC->_colorBurstPhase[1] + pDecodeNTSC->_colorBurstPhase[2] - pDecodeNTSC->_colorBurstPhase[3]);
	bool colorMode = (colorBurstI0*colorBurstI0 + colorBurstQ0*colorBurstQ0) > 2.8 && hf < 16.0f;
	if (colorMode)
		for (int i = 0; i < 4; ++i)
			pDecodeNTSC->_lockedColorBurstPhase[i] = pDecodeNTSC->_colorBurstPhase[i];
	// Color killer hysteresis: We only switch between colour mode and
	// monochrome mode if we stay in the new mode for 128 consecutive
	// lines.
	if (pDecodeNTSC->_colorMode != colorMode) {
		pDecodeNTSC->_hysteresisCount++;
		if (pDecodeNTSC->_hysteresisCount == 128) {
			pDecodeNTSC->_colorMode = colorMode;
			pDecodeNTSC->_hysteresisCount = 0;
		}
	}
	else
		pDecodeNTSC->_hysteresisCount = 0;

	if (pDecodeNTSC->_foundVerticalSync && pDecodeNTSC->_line >= pDecodeNTSC->_topLine && pDecodeNTSC->_line < pDecodeNTSC->_bottomLine) {
		int y0 = pDecodeNTSC->_line - pDecodeNTSC->_topLine;

		// Lines with high amounts of brightness cause more load on the
		// horizontal oscillator which decreases horizontal deflection,
		// causing "blooming" (increase in width).
		int totalSignal = 0;
		for (int i = 0; i < pDecodeNTSC->_active; ++i)
			totalSignal += (int)(readerGetItem(offset + i)) - 60;
		pDecodeNTSC->_crtLoad = 0.4f*pDecodeNTSC->_crtLoad + 0.6f*(pDecodeNTSC->_baseLoad + ((float)totalSignal - 42000.0f)/140000.0f);
		float bloom = clamp(-2.0f, pDecodeNTSC->_bloomFactor*pDecodeNTSC->_crtLoad, 10.0f);
		float horizontalSize = (1.0f - 6.3f*bloom/pDecodeNTSC->_active)*pDecodeNTSC->_horizontalSize;
		float samplesVisible = pDecodeNTSC->_active*horizontalSize;
		float sampleLeft = pDecodeNTSC->_preActive + pDecodeNTSC->_active*(0.5f + pDecodeNTSC->_horizontalPosition - horizontalSize/2.0f);
		pDecodeNTSC->_lefts[y0] = sampleLeft;
		pDecodeNTSC->_widths[y0] = samplesVisible;

		int start = max((int)(sampleLeft) - 10, 0);
		int end = min((int)(sampleLeft + samplesVisible) + 10, pDecodeNTSC->_maxSamplesPerLine - offset);
		int brightness = (int)(pDecodeNTSC->_brightness*100.0 - 7.5f*256.0f*pDecodeNTSC->_contrast)<<8;
		unsigned int* destination = (unsigned int*)(pDecodeNTSC->_dataData + y0*pDecodeNTSC->_dataPitch) + start;

		if (pDecodeNTSC->_colorMode) {
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
			int* p = &pDecodeNTSC->_delay[pDecodeNTSC->_maxSamplesPerLine];
			for (int x = 0; x < 19; ++x)
				p[x] = 0;
			int sp = offset + start;
			for (int x = start; x < end; ++x, --p) {
				// We use a low-pass Finite Impulse Response filter to
				// remove high frequencies (including the color carrier
				// frequency) from the signal. We could just keep a
				// 4-sample running average but that leads to sharp edges
				// in the resulting image.
				int s = (int)(readerGetItem(sp++)) - 60;
				p[0] = s;
				int y = (p[6] + p[0] + ((p[5] + p[1])<<2) + 7*(p[4] + p[2]) + (p[3]<<3))*yContrast + brightness;
				p[6] = s*iqMultipliers[x&3];
				int i = p[12] + p[6] + ((p[11] + p[7])<<2) + 7*(p[10] + p[8]) + (p[9]<<3);
				p[12] = s*iqMultipliers[(x + 3)&3];
				int q = p[18] + p[12] + ((p[17] + p[13])<<2) + 7*(p[16] + p[14]) + (p[15]<<3);
				int r = pDecodeNTSC->_gamma0[clamp(0, (y + 243*i + 160*q)>>16, 255)];
				int g = pDecodeNTSC->_gamma0[clamp(0, (y -  71*i - 164*q)>>16, 255)];
				int b = pDecodeNTSC->_gamma0[clamp(0, (y - 283*i + 443*q)>>16, 255)];
				*(destination++) = 0xff000000 | (r<<16) | (g<<8) | b;
			}
		}
		else {
			int sp = offset + start;
			int yContrast = (int)(pDecodeNTSC->_contrast*46816.0f);
			for (int x = start; x < end; ++x) {
				int s = (int)(readerGetItem(sp++)) - 60;
				int y = pDecodeNTSC->_gamma0[clamp(0, (s*yContrast + brightness)>>16, 255)];
				*(destination++) = 0xFF000000 | (y<<16) | (y<<8) | y;
			}
		}
	}
	offset += pDecodeNTSC->_minSamplesPerLine;
	pDecodeNTSC->_phase = (pDecodeNTSC->_phase + offset)&3;
	consumerRead(offset);

	++pDecodeNTSC->_line;
	if (pDecodeNTSC->_foundVerticalSync && pDecodeNTSC->_line == pDecodeNTSC->_minLinesPerField)
		decodeNTSCpostField(pDecodeNTSC);
}

static void decodeNTSCinit(DecodeNTSC* const pDecodeNTSC)
{
	pDecodeNTSC->_phase=0;
	pDecodeNTSC->_foundVerticalSync=false;
	pDecodeNTSC->_line=0;
	pDecodeNTSC->_baseLoad=0.5;
	pDecodeNTSC->_verticalSync=0;
	pDecodeNTSC->_hysteresisCount=0;
	pDecodeNTSC->_colorMode=false;

	// TODO: make these user-settable
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

	float samplesPerSecond = 157500000.0f/11.0f;
	float us = samplesPerSecond/1000000.0f;  // samples per microsecond

	// Horizontal times in samples.
	float sync = 4.7f*us;
	float breezeway = 0.6f*us;
	pDecodeNTSC->_colorBurstStart = (int)(sync + breezeway);
	float colorBurst = 2.5f*us;
	float backPorch = 1.6f*us;
	float frontPorch = 1.5f*us;
	float blanking = sync + breezeway + colorBurst + backPorch + frontPorch;
	float line = 910.0f;
	pDecodeNTSC->_active = line - blanking;
	pDecodeNTSC->_preActive = blanking - frontPorch;
	// The following parameter specifies how many samples early or late the
	// horizontal sync pulse can be and still be recognized (assuming good
	// signal fidelity). This sets the angle of the diagonal lines that
	// vertical lines become when horizontal sync is lost.
	pDecodeNTSC->_driftSamples = 8;
	pDecodeNTSC->_minSamplesPerLine = (int)(line - (float)pDecodeNTSC->_driftSamples);
	pDecodeNTSC->_maxSamplesPerLine = (int)(line + (float)pDecodeNTSC->_driftSamples);
	// We always consume a scanline at a time, and we won't be called to
	// process until we're sure we have a scanline.
	pDecodeNTSC->_n = pDecodeNTSC->_maxSamplesPerLine;
	pDecodeNTSC->_linePeriod = (int)(line);

	// Vertical times in lines.
	float preSyncLines = 3.0f;
	float syncLines = 3.0f;
	float postSyncLines = 14.0f;
	float lines = 262.5f;
	float blankingLines = preSyncLines + syncLines + postSyncLines;
	float activeLines = lines - blankingLines;
	pDecodeNTSC->_linesVisible = activeLines*pDecodeNTSC->_verticalSize;
	pDecodeNTSC->_lineTop = postSyncLines + activeLines*(0.5f + pDecodeNTSC->_verticalPosition - pDecodeNTSC->_verticalSize/2.0f);
	// The following parameter specifies how many lines early or late the
	// vertical sync pulse can be and still be recognized (assuming good
	// signal fidelity). This sets the "roll speed" of the picture when
	// vertical sync is lost. Empirically determined from video of an IBM
	// 5153 monitor.
	pDecodeNTSC->_driftLines = 14;
	pDecodeNTSC->_minLinesPerField = (int)(lines - (float)pDecodeNTSC->_driftLines);
	pDecodeNTSC->_maxLinesPerField = (int)(lines + (float)pDecodeNTSC->_driftLines);

	for (int i = 0; i < 4; ++i)
		pDecodeNTSC->_colorBurstPhase[i] = pDecodeNTSC->_lockedColorBurstPhase[i] = 0;

	pDecodeNTSC->_crtLoad = pDecodeNTSC->_baseLoad;

	for (int i = 0; i < 256; ++i)
		pDecodeNTSC->_gamma[i] = (int)(pow((float)(i)/255.0f, 1.9f)*255.0f);
	pDecodeNTSC->_gamma0 = &pDecodeNTSC->_gamma[0];
}


static DecodeNTSC s_decodeNTSC;

extern "C" void crtSimInit(unsigned int* pVideoMemoryBGRA)
{
	decodeNTSCinit(&s_decodeNTSC);
	s_decodeNTSC._dataPitch = sizeof(int) * NTSC_SAMPLES_PER_LINE * 2;
	s_decodeNTSC._dataData = (unsigned char*)s_internalTexture;
	s_outputTexture = pVideoMemoryBGRA;
	s_sampleAddIndex = 0;
	s_sampleReadIndex = 0;
}

extern "C" void crtSimAddSample(const unsigned char sample)
{
	s_sampleData[s_sampleAddIndex] = sample;
	s_sampleAddIndex++;
	if (s_sampleAddIndex >= SAMPLEDATA_MAX_SIZE)
	{
		s_sampleAddIndex = 0;
	}
}

extern "C" void crtSimTick(void)
{
	do
	{	
  	decodeNTSCprocess(&s_decodeNTSC);
	} while (s_decodeNTSC._line != 0 || s_decodeNTSC._foundVerticalSync);
}

