#include "AL\al.h"
#include "AL\alc.h"
#include <stdio.h>
#include <string>
#include <math.h>

#ifndef M_PI
#define M_PI atan(1.0)*4
#endif

#define CHANNELS 2
#define SAMPLE_RATE 44100
#define BITS 16
#define LENGTH_S 1
#define T_FRAMES 882
#define NUM_BUFFERS 10

using namespace std;

//This will automatically assume 16-bit audio. It will not work otherwise, will fix this later.
void fillSine(ALvoid *bufferData, short channels, unsigned long frames, int sampleRate, float frequency, short amplitude)
{
	double time = 0;
	short *data = (short*)bufferData, val;

	for(unsigned long i = 0; i < frames; i++)
	{
		time += 1 /(float)sampleRate;
		val = (short)(amplitude * sin(2 * M_PI * frequency * time));

		for(short j = 0; j < channels; j++)
		{
			*data++ = val;
		}
	}
}

int main()
{
	ALCenum error;
	ALCdevice *device;
	device = alcOpenDevice(NULL);
	if(!device)
		printf("Device not found!\n");

	ALboolean enumeration;
	enumeration = alcIsExtensionPresent(NULL, "ALC_ENUMERATION_EXT");
	if (enumeration == AL_FALSE)
		printf("no enumeration supported!\n");

	ALCcontext *context;
	context = alcCreateContext(device, NULL);
	if (!alcMakeContextCurrent(context))
		printf("Context creation failed!");

	ALfloat listenerOri[] = { 0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f };
	alListener3f(AL_POSITION, 0, 0, 1.0f);
	alListener3f(AL_VELOCITY, 0, 0, 0);
	alListenerfv(AL_ORIENTATION, listenerOri);

	ALuint source;
	alGenSources(1, &source);
	alSourcef(source, AL_PITCH, 1);
	alSourcef(source, AL_GAIN, 1);
	alSource3f(source, AL_POSITION, 0, 0, 0);
	alSource3f(source, AL_VELOCITY, 0, 0, 0);
	alSourcei(source, AL_LOOPING, AL_TRUE);

	ALuint buffer;
	alGenBuffers(1, &buffer);

	ALsizei dataSize = T_FRAMES * CHANNELS * (BITS/8);
	ALvoid *bufferData = (ALvoid *)malloc(dataSize);
	fillSine(bufferData, CHANNELS, T_FRAMES, SAMPLE_RATE, 440, 25000);

	//for(int i = 0; i < NUM_BUFFERS; i++)
	//	alBufferData(buffer[i], AL_FORMAT_STEREO16, bufferData, dataSize, SAMPLE_RATE);	
	//alSourceQueueBuffers(source, NUM_BUFFERS, buffer);
	alBufferData(buffer, AL_FORMAT_STEREO16, bufferData, dataSize, SAMPLE_RATE);
	alSourceQueueBuffers(source, AL_BUFFER, &buffer);
	alSourcePlay(source);

	while(true);
	/*short state = 0;
	while(true)
	{
		ALint processed;
		alGetSourcei(source, AL_BUFFERS_PROCESSED, &processed);
		while(processed--)
		{
			alSourceUnqueueBuffers(source, 1, &buffer[state]);
			alSourceQueueBuffers(source, 1, &buffer[state]);

			state++;
			if(state == NUM_BUFFERS)
				state = 0;
		}
	}*/

	// cleanup context
	alDeleteSources(1, &source);
	alDeleteBuffers(1, &buffer);
	device = alcGetContextsDevice(context);
	alcMakeContextCurrent(NULL);
	alcDestroyContext(context);
	alcCloseDevice(device);

	return 0;
}