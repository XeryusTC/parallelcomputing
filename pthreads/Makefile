all: ocr.c
	gcc -Wall -O3 -o ocr ocr.c stringbuffer.c -lm

ex2: ocr-ex2.c stringbuffer.c
	gcc -Wall -O3 -o ocr ocr-ex2.c stringbuffer.c -lm -lpthread

ex3: ocr-ex3.c stringbuffer.c
	gcc -Wall -O3 -o ocr -lm -lpthread ocr-ex3.c stringbuffer.c queue.c
