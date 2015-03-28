all:
	gcc -o fft_test -O2 -ftree-vectorize -I /usr/include/eigen3/ fft_test.cpp -lstdc++

test:
	gcc -o fft_tb -O3  -march=nehalem fft_tb.cpp -lstdc++
