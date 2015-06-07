#ifndef STRING_BUFFER_H
#define STRING_BUFFER_H

#include <pthread.h>

/* Data type for storing 2D greyscale image */
typedef struct imagestruct
{
    int width, height;
    int **imdata;
} *Image;

Image makeImage(int w, int h);
void freeImage(Image im);

typedef struct string {
    char *str;
    int size;
    int end;
} String;

String* newEmptyString(int size);
void freeString(String *str);
void resizeString(String *str);
void appendString(String *str, char c);

typedef struct imagestring {
	char *str;
	Image *im;
	int size, end, front;
	pthread_mutex_t m;
} *ImageString;

ImageString newImageString(int size);
void freeImageString(ImageString imstr);
void resizeImageString(ImageString imstr);
ImageString copyImageString(ImageString org);
void imageStringAppendChar(ImageString imstr, char c);
void imageStringAppendImage(ImageString imstr, Image im);
int imageStringNextPosition(ImageString imstr);

#endif /* STRING_BUFFER_H */
