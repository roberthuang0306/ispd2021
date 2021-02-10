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
