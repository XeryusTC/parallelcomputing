#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    #pragma omp parallel
    {
        printf("Hello ");
        printf("world\n");
    }
    #pragma omp parallel
    printf("Have a nice day!\n");
    printf("Have fun\n");

    return 0;
}
