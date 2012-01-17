#include <malloc.h>
#include <memory.h>

#include <GL/glfw3.h>
#include <GL/glext.h>

const unsigned int WIDTH=910;
const unsigned int HEIGHT=262*2;
const unsigned int MAIN_WINDOW=0;

#define MAX_WINDOWS 8
static GLFWwindow windows[MAX_WINDOWS];
static unsigned char *videoMemory[MAX_WINDOWS];
static GLuint videoTexture[MAX_WINDOWS];

void ShowScreen(unsigned int windowNum, unsigned int w, unsigned int h)
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

void setupGL(unsigned int windowNum, unsigned int w, unsigned int h) 
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

void restoreGL(unsigned int windowNum, unsigned int w, unsigned int h) 
{
	/* Tell OpenGL how to convert from coordinates to pixel values */
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glClearColor(1.0f, 0.f, 1.0f, 1.0f);
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity(); 

	glDisable(GL_TEXTURE_2D);
	glEnable(GL_TEXTURE_RECTANGLE_NV);
	glDisable(GL_DEPTH_TEST);
	(void)windowNum;
}

struct KeyArray
{
	unsigned char lastState;
	unsigned char curState;
	unsigned char debounced;
	GLFWwindow    window;
};

struct KeyArray keyArray[512];

int KeyDown(int key)
{
	return keyArray[key].curState==GLFW_PRESS;
}

int CheckKey(int key)
{
	return keyArray[key].debounced;
}

int CheckKeyWindow(int key, GLFWwindow window)
{
	return keyArray[key].debounced && (keyArray[key].window==window);
}

void ClearKey(int key)
{
	keyArray[key].debounced=0;
}

void kbHandler( GLFWwindow window, int key, int action )
{
	keyArray[key].lastState = keyArray[key].curState;
	keyArray[key].curState = (unsigned char)action;
	keyArray[key].debounced |= (unsigned char)((keyArray[key].lastState==GLFW_RELEASE)&&(keyArray[key].curState==GLFW_PRESS));
	keyArray[key].window = window;
}

void sizeHandler(GLFWwindow window, int xs, int ys)
{
	glfwMakeContextCurrent(window);
	(void)window;
	glViewport(0, 0, xs, ys);
}

void AttachImage(const char* fileName);
void SaveTAP(const char* filename);

uint8_t CheckKeys(int key0, int key1, int key2, int key3, int key4, int key5, int key6, int key7)
{
	uint8_t keyVal=0xFF;

	if (KeyDown(key0))
	{
		keyVal = (uint8_t)(keyVal & ~0x01);
	}
	if (KeyDown(key1))
	{
		keyVal = (uint8_t)(keyVal & ~0x02);
	}
	if (KeyDown(key2))
	{
		keyVal = (uint8_t)(keyVal & ~0x04);
	}
	if (KeyDown(key3))
	{
		keyVal = (uint8_t)(keyVal & ~0x08);
	}
	if (KeyDown(key4))
	{
		keyVal = (uint8_t)(keyVal & ~0x10);
	}
	if (KeyDown(key5))
	{
		keyVal = (uint8_t)(keyVal & ~0x20);
	}
	if (KeyDown(key6))
	{
		keyVal = (uint8_t)(keyVal & ~0x40);
	}
	if (KeyDown(key7))
	{
		keyVal = (uint8_t)(keyVal & ~0x80);
	}

	return keyVal;
}

int windowSetup(void)
{
	unsigned int w;
	unsigned int h;

	/* Initialize GLFW */
	glfwInit(); 

	/* Open screen OpenGL window */
	if(!(windows[MAIN_WINDOW] = glfwOpenWindow((int)WIDTH, (int)HEIGHT, GLFW_WINDOWED, "NTSC Decode", NULL)) ) 
	{ 
		glfwTerminate(); 
		return 1; 
	} 

	glfwSetWindowPos(windows[MAIN_WINDOW], 0, 300);

	glfwMakeContextCurrent(windows[MAIN_WINDOW]);
	setupGL(MAIN_WINDOW, WIDTH, HEIGHT);

	glfwSwapInterval(0);			/* Disable VSYNC */

	glfwGetWindowSize(windows[MAIN_WINDOW], (int*)&w, (int*)&h);

	/* printf("width : %d (%d) , height : %d (%d)\n", w, WIDTH, h, HEIGHT); */
	glfwSetKeyCallback(kbHandler);
	glfwSetWindowSizeCallback(sizeHandler);

	return 0;
}

void windowMainLoop(void)
{
	glfwMakeContextCurrent(windows[MAIN_WINDOW]);
	ShowScreen(MAIN_WINDOW, WIDTH, HEIGHT);
	glfwSwapBuffers();

	glfwPollEvents();
}



