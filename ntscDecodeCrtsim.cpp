#include <stdio.h>
#include <math.h>
#include <memory.h>

typedef unsigned char Sample;  // Sync = 4, Blank = 60, Black = 70, White = 200

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

#define SAMPLE_MAX_SIZE (910*525)
static Sample s_sampleData[SAMPLE_MAX_SIZE];
static int s_sampleOffset = 0;
static int s_sampleAddIndex = 0;
static int s_fieldOffset = 0;

static Sample readerGetItem(const Sample* const pSamples, const unsigned int offset)
{
	unsigned int maxSize = 910+8;
	unsigned int sampleIndex = offset;
	if (sampleIndex >= maxSize)
	{
		sampleIndex = (maxSize-1);
	}
	return pSamples[sampleIndex];
}

static void consumerRead(const int offset)
{
	s_sampleOffset += offset;
	while (s_sampleOffset >= SAMPLE_MAX_SIZE)
	{
		s_sampleOffset -= SAMPLE_MAX_SIZE;
	}
}

struct CompositeMonitor
{
    CompositeMonitor()
		{
        _phase=0;
        _foundVerticalSync=false;
        _line=0;
        _baseLoad=0.5;
        _verticalSync=0;
        _hysteresisCount=0;
        _colorMode=false;

        // TODO: make these user-settable
        _brightness=0.06f;
        _contrast=3.0f;
        _saturation=0.7f;
        _tint=18.0f;
        _horizontalSize=0.95f;
        _horizontalPosition=0;
        _verticalSize=0.93f;
        _verticalPosition=-0.01f;
        _verticalHold=280;
        _horizontalHold=25;
        _bloomFactor=10.0f;

        float samplesPerSecond = 157500000.0f/11.0f;
        float us = samplesPerSecond/1000000.0f;  // samples per microsecond

        // Horizontal times in samples.
        float sync = 4.7f*us;
        float breezeway = 0.6f*us;
        _colorBurstStart = (int)(sync + breezeway);
        float colorBurst = 2.5f*us;
        float backPorch = 1.6f*us;
        float frontPorch = 1.5f*us;
        float blanking = sync + breezeway + colorBurst + backPorch + frontPorch;
        float line = 910.0f;
        _active = line - blanking;
        _preActive = blanking - frontPorch;
        // The following parameter specifies how many samples early or late the
        // horizontal sync pulse can be and still be recognized (assuming good
        // signal fidelity). This sets the angle of the diagonal lines that
        // vertical lines become when horizontal sync is lost.
        _driftSamples = 8;
        _minSamplesPerLine = (int)(line - (float)_driftSamples);
        _maxSamplesPerLine = (int)(line + (float)_driftSamples);
        // We always consume a scanline at a time, and we won't be called to
        // process until we're sure we have a scanline.
        _n = _maxSamplesPerLine;
        _linePeriod = (int)(line);

        // Vertical times in lines.
        float preSyncLines = 3.0f;
        float syncLines = 3.0f;
        float postSyncLines = 14.0f;
        float lines = 262.5f;
        float blankingLines = preSyncLines + syncLines + postSyncLines;
        float activeLines = lines - blankingLines;
        _linesVisible = activeLines*_verticalSize;
        _lineTop = postSyncLines + activeLines*(0.5f + _verticalPosition - _verticalSize/2.0f);
        // The following parameter specifies how many lines early or late the
        // vertical sync pulse can be and still be recognized (assuming good
        // signal fidelity). This sets the "roll speed" of the picture when
        // vertical sync is lost. Empirically determined from video of an IBM
        // 5153 monitor.
        _driftLines = 14;
        _minLinesPerField = (int)(lines - (float)_driftLines);
        _maxLinesPerField = (int)(lines + (float)_driftLines);

        for (int i = 0; i < 4; ++i)
            _colorBurstPhase[i] = _lockedColorBurstPhase[i] = 0;

        _crtLoad = _baseLoad;

        for (int i = 0; i < 256; ++i)
            _gamma[i] = (int)(pow((float)(i)/255.0f, 1.9f)*255.0f);
        _gamma0 = &_gamma[0];
    }

    // We always process a single scanline here.
    void process()
    {
				Sample readerData[910+8];
				memcpy(readerData, &(s_sampleData[s_sampleOffset]), (910+8)*sizeof(Sample));

        // Find the horizontal sync position.
        int offset = 0;
        for (int i = 0; i < _driftSamples*2; ++i, ++offset)
            if ((int)(readerGetItem(readerData, offset)) + (int)(readerGetItem(readerData, offset + 1)) < _horizontalHold*2)
                break;

        // We use a phase-locked loop like real hardware does, in order to
        // avoid losing horizontal sync if the pulse is missing for a line or
        // two, and so that we get the correct "wobble" behavior.
        int linePeriod = _maxSamplesPerLine - offset;
        _linePeriod = (2*_linePeriod + linePeriod)/3;
        _linePeriod = clamp(_minSamplesPerLine, _linePeriod, _maxSamplesPerLine);
        offset = _maxSamplesPerLine - _linePeriod;

        // Find the vertical sync position.
        if (!_foundVerticalSync)
            for (int j = 0; j < _maxSamplesPerLine; j += 57) {
                _verticalSync = ((_verticalSync*232)>>8) + (int)(readerGetItem(readerData, j)) - 60;
                if (_verticalSync < -_verticalHold || _line == 2*_driftLines) {
                    // To render interlaced signals correctly, we need to
                    // figure out where the vertical sync pulse happens
                    // relative to the horizontal pulse. This determines the
                    // vertical position of the raster relative to the screen.
                    _verticalSyncPhase = (float)(j)/(float)(_maxSamplesPerLine);
                    // Now we can find out which scanlines are at the top and
                    // bottom of the screen.
                    _topLine = (int)(0.5f + _lineTop + _verticalSyncPhase);
                    _bottomLine = (int)(1.5f + _linesVisible + _lineTop + _verticalSyncPhase);
                    _line = 0;
                    _foundVerticalSync = true;
                    break;
                }
            }

        // Determine the phase and strength of the color signal from the color
        // burst, which starts shortly after the horizontal sync pulse ends.
        // The color burst is 9 cycles long, and we look at the middle 5
        // cycles.
        int p0 = offset&~3;
        for (int i = _colorBurstStart + 8; i < _colorBurstStart + 28; ++i) {
            static const float colorBurstFadeConstant = 1.0f/128.0f;
            _colorBurstPhase[(i + _phase)&3] =
                _colorBurstPhase[(i + _phase)&3]*(1.0f - colorBurstFadeConstant) +
                (float)((int)(readerGetItem(readerData, p0 + i)) - 60)*colorBurstFadeConstant;
        }
        float total = 0.1f;
        for (int i = 0; i < 4; ++i)
            total += _colorBurstPhase[i]*_colorBurstPhase[i];
        float colorBurstGain = 32.0f/sqrtf(total);
        int phaseCorrelation = (offset + _phase)&3;
        float colorBurstI0 = colorBurstGain*(_colorBurstPhase[2] - _colorBurstPhase[0])/16.0f;
        float colorBurstQ0 = colorBurstGain*(_colorBurstPhase[3] - _colorBurstPhase[1])/16.0f;
        float hf = colorBurstGain*(_colorBurstPhase[0] - _colorBurstPhase[1] + _colorBurstPhase[2] - _colorBurstPhase[3]);
        bool colorMode = (colorBurstI0*colorBurstI0 + colorBurstQ0*colorBurstQ0) > 2.8 && hf < 16.0f;
        if (colorMode)
            for (int i = 0; i < 4; ++i)
                _lockedColorBurstPhase[i] = _colorBurstPhase[i];
        // Color killer hysteresis: We only switch between colour mode and
        // monochrome mode if we stay in the new mode for 128 consecutive
        // lines.
        if (_colorMode != colorMode) {
            _hysteresisCount++;
            if (_hysteresisCount == 128) {
                _colorMode = colorMode;
                _hysteresisCount = 0;
            }
        }
        else
            _hysteresisCount = 0;

        if (_foundVerticalSync && _line >= _topLine && _line < _bottomLine) {
            int y0 = _line - _topLine;
						y0 += s_fieldOffset;

            // Lines with high amounts of brightness cause more load on the
            // horizontal oscillator which decreases horizontal deflection,
            // causing "blooming" (increase in width).
            int totalSignal = 0;
            for (int i = 0; i < _active; ++i)
                totalSignal += (int)(readerGetItem(readerData, offset + i)) - 60;
            _crtLoad = 0.4f*_crtLoad + 0.6f*(_baseLoad + ((float)totalSignal - 42000.0f)/140000.0f);
            float bloom = clamp(-2.0f, _bloomFactor*_crtLoad, 10.0f);
            float horizontalSize = (1.0f - 6.3f*bloom/_active)*_horizontalSize;
            float samplesVisible = _active*horizontalSize;
            float sampleLeft = _preActive + _active*(0.5f + _horizontalPosition - horizontalSize/2.0f);

            int start = max((int)(sampleLeft) - 10, 0);
            int end = min((int)(sampleLeft + samplesVisible) + 10, _maxSamplesPerLine - offset);
            int brightness = (int)(_brightness*100.0 - 7.5f*256.0f*_contrast)<<8;
            unsigned int* destination = (unsigned int*)(_dataData + y0*_dataPitch*2) + start;

            if (_colorMode) {
                int yContrast = (int)(_contrast*1463.0f);
                float radians = (float)(M_PI)/180.0f;
                float tintI = -cosf((103.0f + _tint)*radians);
                float tintQ = sinf((103.0f + _tint)*radians);
                int iqMultipliers[4];
                float colorBurstI = _lockedColorBurstPhase[(2 + phaseCorrelation)&3] - _lockedColorBurstPhase[(0 + phaseCorrelation)&3];
                float colorBurstQ = _lockedColorBurstPhase[(3 + phaseCorrelation)&3] - _lockedColorBurstPhase[(1 + phaseCorrelation)&3];
                iqMultipliers[0] = (int)((colorBurstI*tintI - colorBurstQ*tintQ)*_saturation*_contrast*colorBurstGain*0.352f);
                iqMultipliers[1] = (int)((colorBurstQ*tintI + colorBurstI*tintQ)*_saturation*_contrast*colorBurstGain*0.352f);
                iqMultipliers[2] = -iqMultipliers[0];
                iqMultipliers[3] = -iqMultipliers[1];
                int* p = &_delay[_maxSamplesPerLine];
                for (int x = 0; x < 19; ++x)
                    p[x] = 0;
                int sp = offset + start;
                for (int x = start; x < end; ++x, --p) {
                    // We use a low-pass Finite Impulse Response filter to
                    // remove high frequencies (including the color carrier
                    // frequency) from the signal. We could just keep a
                    // 4-sample running average but that leads to sharp edges
                    // in the resulting image.
                    int s = (int)(readerGetItem(readerData, sp++)) - 60;
                    p[0] = s;
                    int y = (p[6] + p[0] + ((p[5] + p[1])<<2) + 7*(p[4] + p[2]) + (p[3]<<3))*yContrast + brightness;
                    p[6] = s*iqMultipliers[x&3];
                    int i = p[12] + p[6] + ((p[11] + p[7])<<2) + 7*(p[10] + p[8]) + (p[9]<<3);
                    p[12] = s*iqMultipliers[(x + 3)&3];
                    int q = p[18] + p[12] + ((p[17] + p[13])<<2) + 7*(p[16] + p[14]) + (p[15]<<3);
                    int r = _gamma0[clamp(0, (y + 243*i + 160*q)>>16, 255)];
                    int g = _gamma0[clamp(0, (y -  71*i - 164*q)>>16, 255)];
                    int b = _gamma0[clamp(0, (y - 283*i + 443*q)>>16, 255)];
                    *(destination++) = 0xff000000 | (r<<16) | (g<<8) | b;
                }
            }
            else {
                int sp = offset + start;
                int yContrast = (int)(_contrast*46816.0f);
                for (int x = start; x < end; ++x) {
                    int s = (int)(readerGetItem(readerData, sp++)) - 60;
                    int y = _gamma0[clamp(0, (s*yContrast + brightness)>>16, 255)];
                    *(destination++) = 0xFF000000 | (y<<16) | (y<<8) | y;
                }
            }
        }
        offset += _minSamplesPerLine;
        _phase = (_phase + offset)&3;
				consumerRead(offset);

        ++_line;
        if (_foundVerticalSync && _line == _minLinesPerField)
            postField();
    }

    void postField()
    {
        _line = 0;
        _foundVerticalSync = false;
        _crtLoad = _baseLoad;
        _verticalSync = 0;
				s_sampleAddIndex = 0;
				s_sampleOffset = 0;
				s_fieldOffset ^= 1;
    }

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

    int _delay[19+910+8];
    int _gamma[256];
    int* _gamma0;

    int _topLine;
    int _bottomLine;

    int _linePeriod;

    int _hysteresisCount;
    bool _colorMode;

		int _n;
};

static CompositeMonitor s_compositeMonitor;

extern "C" void crtSimInit(unsigned int* pVideoMemoryBGRA)
{
	s_compositeMonitor._dataPitch = 4 * 910;
	s_compositeMonitor._dataData = (unsigned char*)pVideoMemoryBGRA;
	s_sampleAddIndex = 0;
	s_sampleOffset = 0;
	s_fieldOffset = 0;
	memset(s_sampleData, 200, sizeof(Sample)*SAMPLE_MAX_SIZE);
}

extern "C" void crtSimAddSample(const unsigned char sample)
{
	s_sampleData[s_sampleAddIndex] = sample;
	s_sampleAddIndex++;
	if (s_sampleAddIndex >= SAMPLE_MAX_SIZE)
	{
		s_sampleAddIndex = 0;
	}
}

extern "C" void crtSimTick(void)
{
	do
	{	
  	s_compositeMonitor.process();
	} while (s_compositeMonitor._line != 0 || s_compositeMonitor._foundVerticalSync);
	s_sampleAddIndex = 0;
	s_sampleOffset = 0;
	memset(s_sampleData, 200, sizeof(Sample)*SAMPLE_MAX_SIZE);
}

