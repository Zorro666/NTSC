#include <malloc.h>
#include <memory.h>

#include <GL/glfw3.h>
#include <GL/glext.h>

#define MAX_WINDOWS 8
static GLFWwindow windows[MAX_WINDOWS];
static unsigned char *videoMemory[MAX_WINDOWS];
static GLuint videoTexture[MAX_WINDOWS];

static int FindWindowNum(GLFWwindow glfwWindow)
{
	int i;
	for (i = 0; i < MAX_WINDOWS; i++)
	{
		if (windows[i] == glfwWindow)
		{
			return i;
		}
	}
	return -1;
}

static void ShowScreen(unsigned int windowNum, unsigned int w, unsigned int h)
{
	glBindTexture(GL_TEXTURE_RECTANGLE_NV, (GLuint)videoTexture[windowNum]);

	/* glTexSubImage2D is faster when not using a texture range */
	glTexSubImage2D(GL_TEXTURE_RECTANGLE_NV, 0, 0, 0, (GLsizei)w, (GLsizei)h, GL_BGRA, GL_UNSIGNED_INT_8_8_8_8_REV, videoMemory[windowNum]);
	glBegin(GL_QUADS);

	glTexCoord2f(0.0f, 0.0f);
	glVertex2f(-1.0f, 1.0f);

	glTexCoord2f(0.0f, (float)h);
	glVertex2f(-1.0f, -1.0f);

	glTexCoord2f((float)w, (float)h);
	glVertex2f(1.0f, -1.0f);

	glTexCoord2f((float)w, 0.0f);
	glVertex2f(1.0f, 1.0f);
	glEnd();

	glFlush();
}

static void setupGL(unsigned int windowNum, unsigned int w, unsigned int h) 
{
	unsigned int memSize = (unsigned int)w * (unsigned int)h * sizeof(unsigned int);
	videoTexture[windowNum] = windowNum+1;
	videoMemory[windowNum] = (unsigned char*)malloc(memSize);
	memset(videoMemory[windowNum], 0, memSize);

	/* Tell OpenGL how to convert from coordinates to pixel values */
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glClearColor(1.0f, 0.f, 1.0f, 1.0f);
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity(); 

	glDisable(GL_TEXTURE_2D);
	glEnable(GL_TEXTURE_RECTANGLE_NV);
	glBindTexture(GL_TEXTURE_RECTANGLE_NV, videoTexture[windowNum]);

	/* glTextureRangeAPPLE(GL_TEXTURE_RECTANGLE_EXT, 0, NULL); */
	/*	glTexParameteri(GL_TEXTURE_RECTANGLE_EXT, GL_TEXTURE_STORAGE_HINT_APPLE , GL_STORAGE_CACHED_APPLE); */
	/*	glPixelStorei(GL_UNPACK_CLIENT_STORAGE_APPLE, GL_TRUE); */

	glTexParameteri(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);

	glTexImage2D(GL_TEXTURE_RECTANGLE_NV, 0, GL_RGBA, 
							(GLsizei)w, (GLsizei)h, 0, 
							GL_BGRA, GL_UNSIGNED_INT_8_8_8_8_REV, videoMemory[windowNum]);

	glDisable(GL_DEPTH_TEST);
}

struct KeyArray
{
	unsigned char lastState;
	unsigned char curState;
	unsigned char debounced;
	int window;
};

struct KeyArray keyArray[512];

static void kbHandler(GLFWwindow glfwWindow, int key, int action) 
{
	const int window = FindWindowNum(glfwWindow);
	keyArray[key].lastState = keyArray[key].curState;
	keyArray[key].curState = (unsigned char)action;
	keyArray[key].debounced |= (unsigned char)((keyArray[key].lastState == GLFW_RELEASE)&&(keyArray[key].curState == GLFW_PRESS));
	keyArray[key].window = window;
	if (action == GLFW_PRESS)
	{
		printf("Key 0x%X down\n", key);
	}
}

static void sizeHandler(GLFWwindow window, int xs, int ys)
{
	glfwMakeContextCurrent(window);
	glViewport(0, 0, xs, ys);
}

int windowSetup(const unsigned int window, const unsigned int width, const unsigned int height)
{
	unsigned int w;
	unsigned int h;

	/* Initialize GLFW */
	glfwInit(); 

	/* Open screen OpenGL window */
	if(!(windows[window] = glfwOpenWindow((int)width, (int)height, GLFW_WINDOWED, "NTSC Decode", NULL)) ) 
	{ 
		glfwTerminate(); 
		return 1; 
	} 

	glfwSetWindowPos(windows[window], 0, 300);

	glfwMakeContextCurrent(windows[window]);
	setupGL(window, width, height);

	glfwSwapInterval(0);			/* Disable VSYNC */

	glfwGetWindowSize(windows[window], (int*)&w, (int*)&h);

	/* printf("width : %d (%d) , height : %d (%d)\n", w, width, h, height); */
	glfwSetKeyCallback(kbHandler);
	glfwSetWindowSizeCallback(sizeHandler);

	return 0;
}

void windowUpdate(const unsigned int window, const unsigned int width, const unsigned int height)
{
	glfwMakeContextCurrent(windows[window]);
	ShowScreen(window, width, height);
	glfwSwapBuffers();
}

void windowMainLoop(void)
{
	glfwPollEvents();
}

unsigned char* windowGetVideoMemoryBGRA(const unsigned int window)
{
	return videoMemory[window];
}

int windowKeyDown(const int key)
{
	return keyArray[key].curState == GLFW_PRESS;
}

int windowCheckKey(const int key)
{
	return keyArray[key].debounced;
}

int windowCheckKeyWindow(const int key, const int window)
{
	return keyArray[key].debounced && (keyArray[key].window == window);
}

void windowClearKey(const int key)
{
	keyArray[key].debounced=0;
}


