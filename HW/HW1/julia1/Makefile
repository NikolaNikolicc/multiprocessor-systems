# C compiler
CC = gcc
CC_FLAGS = -fopenmp -O2

julia: 
	$(CC) $(CC_FLAGS) julia.c -o julia -lm 

clean:
	rm -f julia julia*.tga
