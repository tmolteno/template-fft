all:
	gcc -o fft_test -O2 -ftree-vectorize -I /usr/include/eigen3/ fft_test.cpp -lstdc++
