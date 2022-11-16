SRC_DIR := src
SRC_FILES := $(shell find $(SRC_DIR) -type f -name '*.cpp')

BUILD_DIR = build
OBJ_DIR := $(BUILD_DIR)/obj
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

LDFLAGS := -lglut -lGLEW -lGL
CXXFLAGS := -Wall

INCLUDE_DIRS = -I$(GLM_INCLUDE)

.PHONY: clean

sph: $(OBJ_FILES)
	g++ -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@ mkdir -p "$(@D)"
	g++ $(CXXFLAGS) $(INCLUDE_DIRS) -c -o $@ $<

clean:
	rm -rf $(BUILD_DIR)