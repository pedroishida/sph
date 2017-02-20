CFLAGS=-Wall -O3
CPPFLAGS=-Wall -O3 --std=c++11 -pthread

TARGETS=sphstar sphdam

build: $(TARGETS)

$(TARGETS): canvas.cpp canvas.h -lSDL2

.PHONY: clean
clean:
	-rm -f $(TARGETS)
