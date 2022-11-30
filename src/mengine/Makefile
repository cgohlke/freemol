all:
	cd src;make

install:
	install -d $(FREEMOL)/bin
	install src/mengine.exe $(FREEMOL)/bin
	install -d $(FREEMOL)/licenses
	install src/LICENSE.txt $(FREEMOL)/licenses/mengine.txt

clean:
	rm -f src/*.o src/mengine.exe
