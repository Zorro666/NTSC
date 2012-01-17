#include <stdio.h>

#include "window.h"

const unsigned int MAIN_WINDOW = 0;

const unsigned int MAIN_WIDTH = 910;
const unsigned int MAIN_HEIGHT = 262*2;

int main(int argc, char* argv[])
{
	int i;
	windowSetup(MAIN_WINDOW, MAIN_WIDTH, MAIN_HEIGHT);

	for (i = 0; i < argc; i++)
	{
		printf("argv[%d] '%s'\n", i, argv[i]);
	}
	while (1)
	{
		unsigned int x;
		unsigned int y;
		unsigned char* const texture = windowGetVideoMemoryBGRA(MAIN_WINDOW);
		i = 0;
		for (y = 0; y < MAIN_HEIGHT; y++)
		{
			for (x = 0; x < MAIN_WIDTH; x++)
			{
				/* BGRA format */
				texture[i+0] = (unsigned char)((x * (y+10) + 100) & 0xFF);
				texture[i+1] = (unsigned char)(((x+10) * y + 50) & 0xFF);
				texture[i+2] = (unsigned char)(((x+5) * (y-10) + 20) & 0xFF);
				texture[i+3] = 0xFF;
				i += 4;
			}
		}
		windowUpdate(MAIN_WINDOW, MAIN_WIDTH, MAIN_HEIGHT);
		windowMainLoop();

		/* 0x100 = ESC */
		if (windowCheckKey(0x100))
		{
			break;
		}
	}
	return -1;
}
