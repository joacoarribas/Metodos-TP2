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
    ~Matriz();

    void multiplicarEscalar(int i); // Esta no se usa por ahora

    void multiplicarMatrices(Matriz& a, Matriz& b);
    void multiplicarVectorDer(vector<double>& x, vector<double>& y);
    void multiplicarVectorIzq(vector<double>& x, vector<double>& y);
    void multiplicarVectoresDameMatriz(std::vector<double>& a, std::vector<double>& b);

    void mas(Matriz& m);
    void menos(Matriz& m);

    vector<double>& operator [] (int fila);

    static double multiplicarVectoresDameValor(std::vector<double>& a, std::vector<double>& b);
    static void cargarVector(std::vector<double>& x);
    static void cerearVector(std::vector<double>& x);

    void etiquetar(int i, int etiqueta);
    void estimar(int i, int etiqueta);
    void mostrar();
    void mostrar2();
    void trasponer(Matriz& traspuesta);
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
