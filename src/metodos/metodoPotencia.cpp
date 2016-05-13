#include <cmath>
#include <algorithm>
#include "../clases/Matriz.h"

#define EPSILON 1.19e-6f

bool igualdadConTolerancia(double a, double b) {
  if (std::abs(a - b) < EPSILON) {
    return true;
  } else {
    return false;
  }
}

double maxAbs(vector<double>& x) {
  double max = std::abs(x[0]);
  int length = x.size();

  for (size_t i = 1; i < length; ++i) {
    if (std::abs(x[i]) > max)
      max = std::abs(x[i]);
  }

  return max;
}

double metodoPotencia(Matriz& A, vector<double>& x) {
  int length = x.size();
  vector<double> y(length, 0);
  int k = 0;
  double c2;

  A.multiplicarVectorDer(x, y);

  double c1 = maxAbs(y);
  double aux = c1;

  do {
    c1 = aux; // La primera vez esto no surte efecto

    for (size_t i = 0; i < length; ++i) // Normalizo vector
      x[i] = y[i] / c1; 

    Matriz::cerearVector(y); // Reseteo el vector y

    A.multiplicarVectorDer(x, y);

    Matriz::cerearVector(x); // Reseteo el vector y

    c2 = maxAbs(y);
    aux = c2;

    ++k;
  
  } while(!igualdadConTolerancia(c1, c2) && k < 15000); // Experimentar con metodo de corte Y Epsilon para la tolerancia

  if (k >= 15000)
    std::cout << "estoy devolviendo fruta" << std::endl;

  //std::cout << "autovalor: " << std::scientific << c2 << std::endl;

  return c2;

}

/*
int main(int argc,char** argv) {
  Matriz A(3,3);

  A[0][0] =  3;
  A[1][0] = -1;
  A[2][0] =  0;

  A[0][1] = -1;
  A[1][1] =  2;
  A[2][1] = -1;

  A[0][2] =  0;
  A[1][2] = -1;
  A[2][2] =  3;


  vector<double> x(3, 1); // x = (1, 1, 1)
 // std::cout << "------------" << std::endl;

 // for (size_t i = 0; i < 3; ++i) 
 //   std::cout << x[i];
 // std::cout << std::endl;

  double res = metodoPotencia(A, x);

  std::cout << "El resultado es: " << std::fixed << res << std::endl;

  std::cout << "El autovector asociado es: " << std::endl;

  for (size_t i = 0; i < 3; ++i) 
    std::cout << x[i] << " ";
  std::cout << std::endl;
  
  
  std::vector<int> y(10);
  for (size_t i = 0; i < 10; ++i)
    y[i] = 9-i;


  for (size_t i = 0; i < 4; ++i) {
    for (size_t i = 0; i < 10; ++i) 
      std::cout << y[i] << " ";
    std::cout << std::endl;
    std::nth_element(y.begin(), y.begin()+i, y.end());
  }

  for (size_t i = 0; i < 10; ++i) 
    std::cout << y[i] << " ";
  std::cout << std::endl;

  
  std::cout << "El maximo es: " << *std::max_element(x.begin(), x.end()) << std::endl; 
  std::cout << " y su posicion es: " << std::distance(x.begin(), std::max_element(x.begin(), x.end())) << std::endl;

  return 0;
}

*/
/*


double metodoPotencia(Matriz& A, vector<double>& x) {
  vector<double> y;
  int length = x.size();
  int k = 0;

  y = A*x;
  double c1 = maxAbs(y);

  for (size_t i = 0; i < length; ++i)
    std::cout << "[ " << std::fixed << y[i] << " ] ";
  std::cout << std::endl;

  for (size_t i = 0; i < length; ++i)
    x[i] = y[i] * (1.0/c1); 

  y = A*x;
  double c2 = maxAbs(y);

  for (size_t i = 0; i < length; ++i)
    std::cout << "[ " << std::fixed << y[i] << " ] ";
  std::cout << std::endl;

  A.mostrar();

  while (!igualdadConTolerancia(c1, c2) && k < 10) {
    std::cout << "c1: " << c1 << std::endl;
    std::cout << "c2: " << c2 << std::endl;

    for (size_t i = 0; i < length; ++i)
      x[i] = y[i] * (1.0/c2); 

    y = A*x;

  for (size_t i = 0; i < length; ++i)
    std::cout << "[ " << std::fixed << y[i] << " ] ";
  std::cout << std::endl;

    double aux = c2;
    c2 = maxAbs(y);
    c1 = aux;

    std::cout << "------------------" << std::endl;
    std::cout << "c1: " << c1 << std::endl;
    std::cout << "c2: " << c2 << std::endl;
    std::cout << "------------------" << std::endl;
    ++k;
  }

  return c2;

}
*/
