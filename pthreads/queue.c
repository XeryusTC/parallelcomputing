#include <assert.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stringbuffer.h"
#include "queue.h"

ImageQueue newImageQueue(int s)
{
	int i;
	ImageQueue iq = malloc(sizeof(struct imagequeue));
	iq->filenames = malloc(sizeof(char*) * s);
	iq->images = malloc(sizeof(Image) * s);
	assert(iq->images != NULL);
	assert(iq->filenames != NULL);
	for (i=0; i<s; ++i)
		iq->filenames[i] = malloc(256 * sizeof(char));
	iq->front = 0;
	iq->end   = 0;
	iq->size  = s;
	return iq;
}

void freeImageQueue(ImageQueue iq)
{
	int i;
	for (i = iq->front; i != iq->end; i = (i + 1) % iq->size)
	{
		free(iq->filenames[i]);
		freeImage(iq->images[i]);
	}
	free(iq->filenames);
	free(iq->images);
	free(iq);
}

int isEmptyImageQueue(ImageQueue iq)
{
	return (iq->front == iq->end);
}

int isFullImageQueue(ImageQueue iq)
{
	return ((iq->end + 1) % iq->size == iq->front);
}

void enqueueImageQueue(ImageQueue iq, Image im, char *str)
{
	strncpy(iq->filenames[iq->end], str, 256);

	iq->images[iq->end] = makeImage(im->width, im->height);
	memcpy(iq->images[iq->end]->imdata[0], im->imdata[0],
			im->width * im->height * sizeof(int));

	iq->end = (iq->end + 1) % iq->size;
}

Image dequeueImageQueue(ImageQueue iq, char *str)
{
	int w, h;
	Image im;

	if (isEmptyImageQueue(iq))
	{
		printf("Image queue is empty\n");
		exit(EXIT_FAILURE);
	}

	w = iq->images[iq->front]->width;
	h = iq->images[iq->front]->height;
	im = makeImage(w, h);
	memcpy(im->imdata[0], iq->images[iq->front]->imdata[0], w*h*sizeof(int));
	freeImage(iq->images[iq->front]);

	strncpy(str, iq->filenames[iq->front], 256);
	iq->front = (iq->front + 1) % iq->size;

	return im;
}
