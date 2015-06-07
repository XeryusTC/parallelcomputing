all: ocr.c
	gcc -Wall -O3 -o ocr ocr.c stringbuffer.c -lm

ex2: ocr.c stringbuffer.c
	gcc -Wall -O3 -o ocr ocr-ex2.c stringbuffer.c -lm -lpthread
