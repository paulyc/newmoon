# newmoon
moon phase calculator from JPL ephems

quick n dirty build system:

make newmoon or just make to build binary

make ephems to download ephemeris data (the 431 is quite large but you don't really need it, just I have the code setup to use it but it should be an option)

make all to do both

then with the binary and the ephems, ./newmoon to calculate your new moons approximately
