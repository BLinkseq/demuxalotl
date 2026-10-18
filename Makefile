VERSION=0.1

all: demuxolotl

demuxolotl:
	@echo "Compiling demuxolotl"
	mkdir -p bin/
	go build -ldflags="-s -w" -o demuxolotl demux.go
