// HelloWorld.c

#include <stdio.h>
#include <omp.h>

#define NUM_THREADS 3

int main(int argc, char *argv[]){   
    int id, nthreads;

    #pragma omp parallel private(id)
    {
        id = omp_get_thread_num();
        nthreads = omp_get_num_threads();
        printf("Zdravo svete! Ja sam %d od %d \n", id, nthreads);  
    }

    return 0;
}

