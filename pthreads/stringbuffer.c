#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
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
