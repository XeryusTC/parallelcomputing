#ifndef QUEUE_H
#define QUEUE_H

#include <pthread.h>
#include "stringbuffer.h"

typedef struct imagequeue {
	Image *images;
	char **filenames;
	int front, end, size;
} *ImageQueue;

ImageQueue newImageQueue(int s);
void freeImageQueue(ImageQueue iq);
int isEmptyImageQueue(ImageQueue iq);
int isFullImageQueue(ImageQueue iq);
void enqueueImageQueue(ImageQueue iq, Image im, char *str);
Image dequeueImageQueue(ImageQueue iq, char *str);

#endif /* QUEUE_H */
