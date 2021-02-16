<<<<<<< HEAD
all:
	mkdir -p out; cd validator; gcc -std=c11 -Wall -O2 ispd_validate.c -o ../out/ispd_validate -lm

=======
GCC= g++
CFLAGS= -g -Wall -std=c++11
EXE= fp
SRC= $(shell find src -type f -name "*.cpp")
OBJ= $(SRC:.cpp=.o)
INCLUDES= -I./src

.PHONY: clean

$(EXE): $(OBJ)
	@echo "Compiling $(EXE)..."
	@$(GCC) $(CFLAGS) $(INCLUDES) -o $(EXE) $(OBJ)

%.o: %.cpp
	@echo "Compiling $<..."
	@$(GCC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	@rm -f $(OBJ) $(EXE) 
>>>>>>> a577815aa3d22cbf40fa63d215d907e48a4e2479
