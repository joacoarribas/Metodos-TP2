CC = g++
CFLAGS = -std=c++11 -o
EFLAGS = -g

BINARIES = main matriz matrizsimetrica gauss descomposicionLU metodoPotencia

.PHONY: main matriz matrizsimetrica gauss descomposicionLU

all: $(BINARIES)

main: main.cpp
	$(CC) main.cpp clases/PLSDA.cpp clases/Matriz.cpp clases/MatrizSimetrica.cpp $(CFLAGS) main

mainTest: mainTest.cpp
	$(CC) mainTest.cpp clases/PLSDA.cpp clases/Matriz.cpp clases/MatrizSimetrica.cpp $(CFLAGS) mainTest

metodoPotencia: metodoPotencia.cpp
	$(CC) metodoPotencia.cpp clases/Matriz.cpp clases/MatrizSimetrica.cpp $(CFLAGS) metodoPotencia

PLSDA:
	$(CC) clases/PLSDA.cpp metodoPotencia.cpp clases/Matriz.cpp clases/MatrizSimetrica.cpp $(CFLAGS) PLSDA

matriz: clases/Matriz.cpp
	$(CC) clases/Matriz.cpp $(CFLAGS) matriz

matrizsimetrica: clases/MatrizSimetrica.cpp
	$(CC) clases/MatrizSimetrica.cpp $(CFLAGS) matrizsimetrica

clean:
	rm -f ./*.o main
	rm -f ./*.o mainTest
	rm -f ./*.o metodoPotencia
	rm -f clases/*.o matriz
	rm -f clases/*.o matrizsimetrica
	rm -f clases/*.o PLSDA
