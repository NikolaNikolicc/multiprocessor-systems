# C compiler
CC = gcc
CC_FLAGS = -fopenmp -O3

julia: 
	$(CC) $(CC_FLAGS) julia.c -o julia -lm 

clean:
	rm -f julia julia*.tga
