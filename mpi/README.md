Since integrate1.c uses sine from math.h, compile it as

`mpicc integrate1.c -lm`

Source - https://stackoverflow.com/questions/16006145/ld-undefined-reference-to-symbol-log2glibc-2-2-5
