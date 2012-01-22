#ifndef NTSC_DECODE_HH
#define NTSC_DECODE_HH

#define NTSC_SAMPLES_PER_LINE (910)
#define NTSC_LINES_PER_FRAME (525)
#define NTSC_LINES_PER_FIELD (262.5)

void ntscDecodeInit(unsigned int* pVideoMemoryBGRA);
void ntscDecodeAddSample(const unsigned char sampleValue);
void ntscDecodeTick(void);

extern int windowCheckKey(const int key);
extern void windowClearKey(const int key);


#endif /* #ifndef NTSC_DECODE_HH */
