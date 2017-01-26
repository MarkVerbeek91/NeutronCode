#!/bin/bash

#make clean
#make
#cd Analitics 
#./bin/Debug/NeutronCode.exe
DATE=$(date +"%Y-%m-%d-%H-%M")
NAME=Analitics/output-$DATE.png
gprof ./bin/Debug/NeutronCode.exe | gprof2dot | dot -Tpng -o $NAME

#cd ..
