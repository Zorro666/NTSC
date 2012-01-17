#include <stdio.h>

#include "window.h"

int main(int argc, char* argv[])
{
	int i;
	windowSetup();

	for (i = 0; i < argc; i++)
	{
		printf("argv[%d] '%s'\n", i, argv[i]);
	}
	while (1)
	{
		windowMainLoop();
	}
	return -1;
}
