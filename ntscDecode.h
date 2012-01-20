#ifndef NTSC_DECODE_HH
#define NTSC_DECODE_HH

#define NTSC_SAMPLES_PER_LINE (910)
#define NTSC_LINES_PER_FIELD (262)
#define NTSC_FIELDS_PER_IMAGE (2)

void ntscDecodeInit(unsigned int* pVideoMemoryBGRA);
void ntscDecodeAddSample(const unsigned char sampleValue);
void ntscDecodeTick(void);

#endif /* #ifndef NTSC_DECODE_HH */
