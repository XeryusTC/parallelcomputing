#ifndef STRING_BUFFER_H
#define STRING_BUFFER_H

typedef struct string {
    char *str;
    int size;
    int end;
} String;

String* newEmptyString(int size);
void freeString(String *str);
void resizeString(String *str);
void appendString(String *str, char c);

#endif /* STRING_BUFFER_H */
