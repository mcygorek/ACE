
all: 
	cd src && $(MAKE) && cd -
	cd src/tools && $(MAKE) && cd -

.PHONY: clean

clean:
	rm bin/* 2>/dev/null || true
	rm tools/* 2>/dev/null || true

