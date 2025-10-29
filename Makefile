# ==== Config ====
CC      := gcc
CFLAGS  := -Wall -Wextra -O2 -I/src
LDLIBS  := -lm -llapacke -lopenblas

SRC_DIR := tests/src
BIN_DIR := tests/bin
LIB_SRC := $(wildcard ./src/*.c)
TEST_SRC := $(wildcard $(SRC_DIR)/*.c)
BIN := $(patsubst $(SRC_DIR)/%.c,$(BIN_DIR)/%,$(TEST_SRC))

.PHONY: all clean run run-one

all: $(BIN)

$(BIN_DIR)/%: $(SRC_DIR)/%.c $(LIB_SRC) | $(BIN_DIR)
	$(CC) $(CFLAGS) $< $(LIB_SRC) -o $@ $(LDLIBS)

$(BIN_DIR):
	@mkdir -p $(BIN_DIR)

tests: $(BIN)
	@set -e; \
	for exe in $(BIN); do \
		echo "=== $$exe ==="; \
		"$$exe"; \
		echo; \
	done

tests-one: $(BIN)
	@test -n "$(NAME)" || (echo "Use: make tests-one NAME=arquivo"; exit 1)
	@echo "=== $(BIN_DIR)/$(NAME) ==="
	@$(BIN_DIR)/$(NAME)

clean:
	rm -rf $(BIN_DIR)
