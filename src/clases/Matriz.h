#ifndef MATRIZ_H_
#define MATRIZ_H_

class MatrizSimetrica;
#include "MatrizSimetrica.h"
#include <vector>
#include <assert.h>
#include <cstdlib>
#include <time.h>
#include <iostream>

using namespace std;

class Matriz {

  public:

    Matriz(int filas, int columnas);

    Matriz& operator * (int i);
    Matriz& operator * (Matriz& m);
    Matriz& operator * (MatrizSimetrica& m);
    vector<double>& operator * (vector<double>& m);

    Matriz& operator + (Matriz& m);
    Matriz& operator - (Matriz& m);

    vector<double>& operator [] (int fila);

    Matriz& transponer();
    void randomizar(int semilla);
    void mostrar();

    int dimensionFilas();
    int dimensionColumnas();

  private:

	/*******************************
	 *          Variables          *
	 *******************************/ 

    int filas, columnas;
    vector< vector<double> > matriz;

	/*******************************
	 *          Funciones          *
	 *******************************/ 

    void clean();
};

#endif
