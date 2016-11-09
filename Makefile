# makefile for version2_thread
# from fileDisposal.o , the command is automatically derived by make

Object =  coordinaryTransform.o    fileDisposal.o    initialization.o \
          main.o    matrixFunction.o    MOPSOAidFunction.o    \
          MOPSOFunction.o    parameterSettings.o

GL : $(Object)
	g++ -o GL $(Object) -lpthread

coordinaryTransform.o : MOPSO.h
	g++ -c coordinaryTransform.cpp

fileDisposal.o : MOPSO.h

initialization.o : MOPSO.h

main.o : MOPSO.h

matrixFunction.o : MOPSO.h

MOPSOAidFunction.o : MOPSO.h

MOPSOFunction.o : MOPSO.h

parameterSettings.o : MOPSO.h

.PHONY : clear
clear :
	-rm $(Object) GL

.PHONY : delete
delete :
	-rm -r /home/ws/zzZyj/MOPSO/data/answer/newAnswer
