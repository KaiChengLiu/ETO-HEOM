.PHONY: all clean

all:
	$(MAKE) -C etoheom-cpu/lib
	$(MAKE) -C etoheom-cpu/src
	$(MAKE) -C etoheom-gpu/lib
	$(MAKE) -C etoheom-gpu/src

clean:
	$(MAKE) -C etoheom-gpu/src clean
	$(MAKE) -C etoheom-gpu/lib clean
	$(MAKE) -C etoheom-cpu/src clean
	$(MAKE) -C etoheom-cpu/lib clean

