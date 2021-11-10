gcc_nautyc = nauty27rc2/nauty.c nauty27rc2/nautil.c nauty27rc2/naugraph.c nauty27rc2/schreier.c nauty27rc2/naurng.c
gcc_ignore_errors = -Wno-write-strings
gcc_flags = -std=c++11 -O4 -march=native -fopenmp -g
gcc_compile = g++ ${gcc_flags} ${gcc_nautyc} ${gcc_ignore_errors}

generating_algorithm: generating_algorithm.cpp
	${gcc_compile} generating_algorithm.cpp -o generating_algorithm
