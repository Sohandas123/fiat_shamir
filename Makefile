#make file - this is a comment section

CC=gcc  #compiler
CFLAGS=-I.
DEPS = fiat_shamir.h
OBJ = fiat_shamir.o functions.o
TARGET=output

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(TARGET):	$(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)
	

clean:
	rm *.o $(TARGET)
