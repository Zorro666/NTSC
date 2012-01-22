#ifndef NTSC_DECODE_HH
#define NTSC_DECODE_HH

#define NTSC_SAMPLES_PER_LINE (910)
#define NTSC_LINES_PER_FRAME (525)
#define NTSC_LINES_PER_FIELD (262.5)

#define DISPLAY_RGB 		(0)
#define DISPLAY_Y 			(1)
#define DISPLAY_CHROMA 	(2)
#define DISPLAY_I 			(3)
#define DISPLAY_Q 			(4)
#define DISPLAY_SIGNAL 	(5)
#define DISPLAY_RED 		(6)
#define DISPLAY_GREEN 	(7)
#define DISPLAY_BLUE 		(8)
#define DISPLAY_MAX 		(8)

void ntscDecodeInit(unsigned int* pVideoMemoryBGRA);
void ntscDecodeAddSample(const unsigned char sampleValue);
void ntscDecodeTick(void);

extern int windowCheckKey(const int key);
extern void windowClearKey(const int key);


#endif /* #ifndef NTSC_DECODE_HH */
