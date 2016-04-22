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
    vector<double>& operator * (vector<double>& m);
    Matriz& operator * (Matriz& m);
    Matriz& operator * (MatrizSimetrica& m);

    Matriz& operator + (Matriz& m);

    Matriz& operator - (Matriz& m);

    vector<float>& operator [] (int fila);

    void randomizar(int semilla);
    void transponer();
    void mostrar();

    int dimensionFilas();
    int dimensionColumnas(); 

  private:

	/*******************************
	 *          Variables          *
	 *******************************/ 

    int filas, columnas;
    vector< vector<float> > matriz;

	/*******************************
	 *          Funciones          *
	 *******************************/ 

    void clean();
};

#endif
