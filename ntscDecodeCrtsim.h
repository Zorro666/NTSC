#ifndef NTSC_DECODE_CRTSIM_HH
#define NTSC_DECODE_CRTSIM_HH

void crtSimInit(unsigned int* pVideoMemoryBGRA);
void crtSimAddSample(const unsigned char sampleValue);
void crtSimTick(const int displayModeFlags);

#endif /* #ifndef NTSC_DECODE_CRTSIM_HH */
