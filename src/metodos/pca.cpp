#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "metodoPotencia.cpp"

void calcular_media(Matriz& matriz, vector<double>& v) { //ver si se pierde precisión con la división
  int n = matriz.dimensionFilas();
  int m = matriz.dimensionColumnas();
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i) {
      v[j] = (v[j] + matriz[i][j]); //acumulo
    }
    v[j] = v[j] / (double)n; //calculo promedio final
  }
  std::cout << "8" << std::endl;
}

void calcular_matriz_covarianza(Matriz& matriz, Matriz& res){
  int n = matriz.dimensionFilas();
  int m = matriz.dimensionColumnas();

  vector<double> media(m);
  
  calcular_media(matriz, media);
  
  std::cout << "2" << std::endl;
  double divisor = sqrt ((double)n - 1.0);

  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i) {
      res[i][j] = (matriz[i][j] - media[i]) / divisor; //escribo resultado en la matriz
    }
  std::cout << "3" << std::endl;
  }
}

void calcular_base_ortonormal(Matriz& m, Matriz& matriz_ortonormal, int alfa){ //deja en matriz_ortonormal una matriz de alfa columnas
  double autovalor;
  int n = matriz_ortonormal.dimensionFilas();
  vector<double> aux(n);  

  for (int i = 0; i < alfa; ++i) { //repito alfa veces (hay que experimentar con dicho valor)
    Matriz::cargarVector(aux);
  std::cout << "por entrar a metodo potencia" << std::endl;

    autovalor = metodoPotencia(m, aux); //calculo el i-ésimo autovalor, en aux queda el autovector
  std::cout << "5" << std::endl;
    for (int k = 0; k < n; ++k){
      matriz_ortonormal[i][k] = aux[k]; //completo matriz con base ortonormal
  std::cout << "6" << std::endl;
    }
    /* Deflación */
    Matriz auxiliar(n, n);
    auxiliar.multiplicarVectoresDameMatriz(aux, aux);
    auxiliar.multiplicarEscalar(autovalor);
    m.menos(auxiliar);
  } 
  /* Tengo en matriz_ortonormal la matriz con base de autovectores de matriz */
}

void PCA(Matriz& matriz, Matriz& res, int alfa){
  int n = matriz.dimensionFilas();
  int m = matriz.dimensionColumnas();

  Matriz X(n, m);
  std::cout << "1" << std::endl;
  calcular_matriz_covarianza(matriz, X);
  std::cout << "4" << std::endl;

  Matriz Xt(m, n);
  X.trasponer(Xt);
  
  Matriz m_covarianza(n, n); //revisar dimensiones

  X.multiplicarMatrices(Xt, m_covarianza); //crear matriz covarianza

  Matriz m_ortonormal(alfa, m);

  /* Reducción de la dimensión */
  calcular_base_ortonormal(m_covarianza, m_ortonormal, alfa); //deja en matriz_ortonormal una matriz de alfa columnas
  std::cout << "7" << std::endl;
  /* Transformación característica */
  Matriz m_ortonormal_traspuesta(alfa, n);
  std::cout << "8" << std::endl;
  m_ortonormal.trasponer(m_ortonormal_traspuesta);
  std::cout << "9" << std::endl;

  /* Transformación característica */
  res.multiplicarMatrices(m_ortonormal_traspuesta, matriz); 
  std::cout << "10" << std::endl;

}


//void calcular_matriz(Matriz& matriz) { NO ME BORREN ESTO POR AHORA PORFIS
//  int n = matriz.dimensionFilas();
//  int m = matriz.dimensionColumnas();
//
//  vector<double> media(n);
//  calcular_media(matriz, media);
//  double divisor = sqrt ((double)n - 1.0);
//  for (int i = 0; i < n; ++i) {
//    for (int j = 0; j < m; ++j) {
//      matriz[i][j] = (matriz[i][j] - media[i]) / divisor;
//    }
//  }
//  /* Cuando termina este paso tengo la matriz de covarianza en matriz */
//  //ahora tengo que diagonalizar matriz
//  Matriz *base_ortonormal = new Matriz(n, n);
//  vector<double> aux(base_ortonormal->dimensionColumnas());  
//  int alfa = 10; //defino valor arbitrario
//  double autovalor;
//
//  for (int i = 0; i < alfa; ++i) {
//    Matriz::cargarVector(aux);
//    autovalor = metodoPotencia(matriz, aux); //calculo el i-ésimo autovalor
//    for (int k = 0; k < base_ortonormal->dimensionColumnas(); ++k){
//      (*base_ortonormal)[i][k] = aux[k]; //completo matriz con base ortonormal
//    }
//    /* Deflación */
//    Matriz auxiliar = Matriz::multiplicarVectoresDameMatriz(aux, aux);
//    matriz- (auxiliar * autovalor);
//  } 
//  /* Tengo en matriz_ortonormal la matriz con base de autovectores de matriz */
//  
//}

