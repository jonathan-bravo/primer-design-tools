CXX      := g++
PRIMER3  ?= ./primer3

CXXFLAGS := -O3 -std=c++17 -I$(PRIMER3)/src -pthread \
            -DPRIMER3_PATH=\"$(PRIMER3)\"

# Detect macOS to handle Apple Clang strictness warnings
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    CXXFLAGS += -Wno-c++11-narrowing
endif

DEPFLAGS  = -MMD -MP

P3_OBJS := \
    $(PRIMER3)/src/libprimer3.o    \
    $(PRIMER3)/src/thal_primer.o   \
    $(PRIMER3)/src/dpal_primer.o   \
    $(PRIMER3)/src/oligotm.o       \
    $(PRIMER3)/src/p3_seq_lib.o    \
    $(PRIMER3)/src/read_boulder.o  \
    $(PRIMER3)/src/print_boulder.o \
    $(PRIMER3)/src/format_output.o \
    $(PRIMER3)/src/masker.o        \
    $(PRIMER3)/src/amplicontm.o

LDLIBS := -lm

BIN_DIR   := bin
BUILD_DIR := build
SRC_DIR   := src

TARGET := $(BIN_DIR)/dpro

SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

all: $(TARGET)

$(TARGET): $(OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ $(P3_OBJS) -o $@ $(LDLIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

-include $(DEPS)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BIN_DIR) $(BUILD_DIR)

.PHONY: all clean