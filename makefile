

.PHONY: all unix test clean xcode

all: unix



xcode:
	#@echo "making... ${JOBS} $(MAKEFLAGS) "
	@if [ ! -d "./build_xcode" ]; then mkdir build_xcode; fi
	@cmake -G Xcode -S ./ -B ./build_xcode || cmake3 -G Xcode -S ./ -B ./build_xcode
	#@make -j4 -C ./build
	#@if [ ! -d "./bin" ]; then mkdir bin; fi
	#@cp -p ./build/bin/*  ./bin
	


unix:  CMakeLists.txt
	#@echo "making... ${JOBS} $(MAKEFLAGS) "
	@if [ ! -d "./build" ]; then mkdir build; fi
	@cmake -S ./ -B ./build -DBUILD_TESTING=OFF || cmake3 -S ./ -B ./build -DBUILD_TESTING=OFF
	@make -j4 -C ./build
	@if [ ! -d "./bin" ]; then mkdir bin; fi
	@cp -p ./build/bin/groot ./bin/groot


test: CMakeLists.txt
	@cmake -S ./ -B ./build -DBUILD_TESTING=ON || cmake3 -S ./ -B ./build -DBUILD_TESTING=ON
	@make -j4 -C ./build
	@ctest --test-dir ./build --output-on-failure


clean: 
	@echo "cleaning..."
	@if [ -d "./build" ]; then rm -rf build; fi
	@if [ -d "./build_xcode" ]; then rm -rf build_xcode; fi
	@if [ -d "./bin" ]; then rm -rf bin; fi
	
