export LDFLAGS="-L/opt/libstdc++-compat -Wl,-rpath,/opt/libstdc++-compat"
cmake -DCMAKE_EXE_LINKER_FLAGS="-L/opt/libstdc++-compat -Wl,-rpath,/opt/libstdc++-compat" ..
make

