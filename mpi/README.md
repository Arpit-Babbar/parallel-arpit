Since integrate1.c uses sine from math.h, compile it as

<<<<<<< HEAD
`mpicc integrate1.c -lm`
=======
`mpicc hello.c -lm`
>>>>>>> 1e2a275... Update README.md

Source - https://stackoverflow.com/questions/16006145/ld-undefined-reference-to-symbol-log2glibc-2-2-5
