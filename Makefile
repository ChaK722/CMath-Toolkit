CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -O3 -Iinclude

SRC_DIR = src
OBJ_DIR = build
TARGET = cmath_test

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRCS))

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@
	@echo "Build complete: execute file ./$(TARGET)"

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR) $(TARGET)
	@echo "Clean complete"

.PHONY: all clean