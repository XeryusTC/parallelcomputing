#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "stringbuffer.h"

String* newEmptyString(int size)
{
    assert(size > 0);
    String *str = malloc(sizeof(struct string));
    assert(str != NULL);
    str->str = malloc(sizeof(char) * size);
    assert(str->str != NULL);
    str->size = size;
    str->end = 0;
    str->str[0] = '\0';
    return str;
}

void freeString(String *str)
{
    free(str->str);
    free(str);
}

void resizeString(String *str)
{
    str->size *= 2;
    str->str = realloc(str->str, str->size * sizeof(char));
    assert(str->str != NULL);
}

void appendString(String *str, char c)
{
    if (str->end == str->size - 2)
        resizeString(str);
    str->str[str->end++] = c;
    str->str[str->end] = '\0';
}

ImageString newImageString(int size)
{
	assert(size > 0);
	ImageString imstr = malloc(sizeof(struct imagestring));
	assert(imstr != NULL);
	imstr->str = malloc(sizeof(char)  * size);
	imstr->im  = malloc(sizeof(Image) * size);
	assert(imstr->str != NULL);
	assert(imstr->im  != NULL);
	imstr->size = size;
	imstr->end = 0;
	imstr->str[0] = '\0';
	imstr->front = 0;
	pthread_mutex_init(&imstr->m, NULL);
	return imstr;
}

void freeImageString(ImageString imstr)
{
	pthread_mutex_destroy(&imstr->m);
	/* clean images if they haven't got an associated character yet */
	free(imstr->im);
	free(imstr->str);
	free(imstr);
}

void resizeImageString(ImageString imstr)
{
	imstr->size *= 2;
	imstr->str = realloc(imstr->str, imstr->size * sizeof(char));
	imstr->im  = realloc(imstr->im,  imstr->size * sizeof(Image));
	assert(imstr->str != NULL);
	assert(imstr->im  != NULL);
}

void imageStringAppendChar(ImageString imstr, char c)
{
	pthread_mutex_lock(&imstr->m);
	if (imstr->end == imstr->size - 2)
		resizeImageString(imstr);
	imstr->str[imstr->end++] = c;
	imstr->str[imstr->end] = '\0';
	pthread_mutex_unlock(&imstr->m);
}

void imageStringAppendImage(ImageString imstr, Image im)
{
	pthread_mutex_lock(&imstr->m);
	if (imstr->end == imstr->size - 2)
		resizeImageString(imstr);
	imstr->im[imstr->end] = makeImage(im->width, im->height);
	memcpy(imstr->im[imstr->end]->imdata[0], im->imdata[0],
			im->width * im->height * sizeof(int));
	imstr->str[imstr->end++] = '_';
	imstr->str[imstr->end] = '\0';
	pthread_mutex_unlock(&imstr->m);
}

int imageStringNextPosition(ImageString imstr)
{
	int pos;
	pthread_mutex_lock(&imstr->m);
	if (imstr->front == imstr->end)
	{
		printf("Tried to retrieve an ImageString position beyond the end of"
			   "the string\n");
		pthread_mutex_unlock(&imstr->m);
		return -1;
	}
	pos = imstr->front++;
	pthread_mutex_unlock(&imstr->m);
	return pos;
}
