#ifndef PLSDA_H_
#define PLSDA_H_

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "Matriz.h"
#include "../metodoPotencia.cpp"

class PLSDA {

  public:
    PLSDA();
    Matriz& transformacionCaracteristica(Matriz m, vector< vector<double> > w, int n, int dimensiones);
    Matriz& PLSDAMethod(Matriz valores, int dimensiones);
    double multiplicarVectores(vector<double> a, vector<double> b);
    void cargarVector(vector<double> &x);
    void normalizar(vector<double> &x);
};

#endif