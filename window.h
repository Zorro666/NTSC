#ifndef WINDOW_HH
#define WINDOW_HH

extern int windowSetup(const unsigned int window, const unsigned int width, const unsigned int height);

extern void windowUpdate(const unsigned int window, const unsigned int width, const unsigned int height);
extern void windowMainLoop(void);

extern unsigned char* windowGetVideoMemoryBGRA(const unsigned int window);

int windowCheckKey(const int key);

#endif /* #ifndef WINDOW_HH */
