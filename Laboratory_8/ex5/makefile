CC=nvcc

LIBS=
INCLUDES=-I../../
LIB_FLAGS=-lm

MAIN_NAME=elementWiseAddition_ex5
MAIN_BIN=$(MAIN_NAME)
MAIN_SRC=$(MAIN_NAME).cu

all: $(MAIN_BIN)

$(MAIN_BIN): $(MAIN_SRC)
	$(CC) $^ -o $(MAIN_BIN) $(LIBS) $(INCLUDES) $(LIB_FLAGS)

clean:
	rm -rf $(MAIN_BIN)_*
