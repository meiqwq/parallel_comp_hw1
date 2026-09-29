# We will benchmark you against Intel MKL implementation, the default processor vendor-tuned implementation.
# This makefile is intended for the Intel C compiler.
# Your code must compile (with icc) with the given CFLAGS. You may experiment with the OPT variable to invoke additional compiler options.

CC = icc
OPT = 
CFLAGS = -Wall -DGETTIMEOFDAY -std=c99 $(OPT) -diag-disable=10441
LDFLAGS = -Wall 
# mkl is needed for blas implementation
LDLIBS = -qmkl=sequential -lpthread -lm -diag-disable=10441

targets = benchmark-final benchmark-naive benchmark-blas benchmark-blocked \
          benchmark-test-final benchmark-test-naive benchmark-test-blas benchmark-test-blocked
objects = benchmark.o benchmark-test.o sgemm-final.o sgemm-naive.o sgemm-blas.o sgemm-blocked.o

.PHONY : default
default : all

.PHONY : all
all : clean $(targets)

benchmark-final : benchmark.o sgemm-final.o
	$(CC) -o $@ $^ $(LDLIBS)
benchmark-naive : benchmark.o sgemm-naive.o
	$(CC) -o $@ $^ $(LDLIBS)
benchmark-blas : benchmark.o sgemm-blas.o
	$(CC) -o $@ $^ $(LDLIBS)
benchmark-blocked : benchmark.o sgemm-blocked.o
	$(CC) -o $@ $^ $(LDLIBS)

benchmark-test-final : benchmark-test.o sgemm-final.o
	$(CC) -o $@ $^ $(LDLIBS)
benchmark-test-naive : benchmark-test.o sgemm-naive.o
	$(CC) -o $@ $^ $(LDLIBS)
benchmark-test-blas : benchmark-test.o sgemm-blas.o
	$(CC) -o $@ $^ $(LDLIBS)
benchmark-test-blocked : benchmark-test.o sgemm-blocked.o
	$(CC) -o $@ $^ $(LDLIBS)

%.o : %.c
	$(CC) -c $(CFLAGS) $<

.PHONY : clean
clean:
	rm -f $(targets) $(objects)
