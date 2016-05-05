#ifndef METODOPOTENCIA_H_
#define METODOPOTENCIA_H_

#include <cmath>
#include <algorithm>
#include "clases/Matriz.h"

class met{
  public:

    static bool igualdadConTolerancia(double a, double b);

    static double maxAbs(vector<double>& x);

    static double metodoPotencia(Matriz& A, vector<double>& x);
};

#endif
