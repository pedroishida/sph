CFLAGS=-Wall -O3
CPPFLAGS=-Wall -O3

TARGETS=sphstar sphdam

build: $(TARGETS)

$(TARGETS): canvas.cpp -lSDL2

.PHONY: clean
clean:
	-rm -f $(TARGETS)
