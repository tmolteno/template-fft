# 
#   Huild the FFT template testbench
# 
#   Author: Tim Molteno, tim@physics.otago.ac.nz
# 
#   Copyright Tim Molteno, 2008-2015.
#   Licensed under the GPL v3.
# 
all:
	gcc -o fft_test -O2 -ftree-vectorize -I /usr/include/eigen3/ fft_test.cpp -lstdc++

test:
	gcc -o fft_tb -O3 -march=native fft_tb.cpp -lstdc++
	./fft_tb
