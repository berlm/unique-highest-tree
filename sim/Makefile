CC ?= gcc
CFLAGS ?= -O2 -Wall -Wextra -std=c99

mapping_sim: mapping_sim.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f mapping_sim

.PHONY: clean
