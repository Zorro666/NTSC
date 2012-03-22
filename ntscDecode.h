#ifndef NTSC_DECODE_HH
#define NTSC_DECODE_HH

#define NTSC_SAMPLES_PER_LINE (910)
#define NTSC_LINES_PER_FRAME (525)
#define NTSC_LINES_PER_FIELD (262.5)

#define DISPLAY_RGB 				(0)
#define DISPLAY_Y 					(1)
#define DISPLAY_CHROMA 			(2)
#define DISPLAY_I 					(3)
#define DISPLAY_Q 					(4)
#define DISPLAY_SIGNAL 			(5)
#define DISPLAY_RED 				(6)
#define DISPLAY_GREEN 			(7)
#define DISPLAY_BLUE 				(8)
#define DISPLAY_MAX 				(9)

#define DISPLAY_INTERLACED	(1 << 0)
#define DISPLAY_LEE_MODE		(1 << 1)

void ntscDecodeInit(unsigned int* pVideoMemoryBGRA);
void ntscDecodeAddSample(const unsigned char sampleValue);
void ntscDecodeTick(void);

void ntscEncodeInit(unsigned char* pOutputSignal, const unsigned int outputSize);
unsigned char ntscEncodeAddSample(const unsigned int RGB);
void ntscEncodeTick(void);

extern int windowCheckKey(const int key);
extern void windowClearKey(const int key);

extern inline float clampFloat(const float value, const float low, const float high) 
{ 
	if (value < low)
	{
		return low;
	}
	else if (value > high)
	{
		return high;
	}
	return value;
}

extern inline int minInt(const int a, const int b)
{
	if (a < b)
	{
		return a;
	}
	return b;
}

extern inline int maxInt(const int a, const int b)
{
	if (a > b)
	{
		return a;
	}
	return b;
}

extern inline int clampInt(const int value, const int low, const int high) 
{ 
	if (value < low)
	{
		return low;
	}
	else if (value > high)
	{
		return high;
	}
	return value;
}

#endif /* #ifndef NTSC_DECODE_HH */
