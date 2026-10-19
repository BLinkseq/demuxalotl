VERSION=0.1

all: demuxolotl

demuxolotl:
	@echo "Compiling demuxolotl"
	mkdir -p bin/
	rm -f bin/demuxolotl
	go build -ldflags="-s -w" -o bin/demuxolotl demuxolotl.go decoder.go io.go pamld.go
