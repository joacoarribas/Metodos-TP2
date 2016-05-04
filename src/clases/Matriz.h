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

    void etiquetar(int i, int etiqueta);
    void estimar(int i, int etiqueta);
    void mostrar();
    Matriz& transponer();
    void randomizar(int semilla);

    int dameEtiqueta(int i);
    int dameEstimacion(int i);
    int dimensionFilas();
    int dimensionColumnas();

  private:

	/*******************************
	 *          Variables          *
	 *******************************/ 

    int filas, columnas;
    vector< vector<double> > matriz;
    vector<int> etiquetas;
    vector<int> estimaciones;

	/*******************************
	 *          Funciones          *
	 *******************************/ 

    void clean();
};

#endif
