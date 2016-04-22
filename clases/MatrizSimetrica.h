#ifndef MATRIZSIMETRICA_H_
#define MATRIZSIMETRICA_H_

class Matriz;
#include "Matriz.h"
#include <vector>
#include <assert.h>
#include <cstdlib>
#include <time.h>
#include <iostream>

using namespace std;

class MatrizSimetrica {

  public:

    MatrizSimetrica(int filas, int columnas);

    MatrizSimetrica& operator * (int i);
    Matriz& operator * (MatrizSimetrica& m);
    Matriz& operator * (Matriz& m);

    void  set(int, int, float);
    float get(int, int);

    void randomizar(int semilla);
    void transponer();
    void mostrar();

    int  dimensionFilas();
    int  dimensionColumnas();
    
    void setTriangular(bool esTriangular);
    bool esTriangular();

  private:

	/*******************************
	 *          Variables          *
	 *******************************/ 

    int filas, columnas;
    vector< vector<float> > matriz;
    bool transpuesta;
    bool triangular;

	/*******************************
	 *          Funciones          *
	 *******************************/ 

    void clean();
};

#endif
