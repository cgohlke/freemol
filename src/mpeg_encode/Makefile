
all:
	cd src;make

install:
	install -d $(FREEMOL)/bin
	install src/mpeg_encode.exe $(FREEMOL)/bin
	install -d $(FREEMOL)/licenses
	install src/COPYRIGHT $(FREEMOL)/licenses/mpeg_encode.txt

clean:
	rm -f src/*.o src/mpeg_encode.exe
