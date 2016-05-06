#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "metodoPotencia.cpp"

void calcular_media(Matriz& matriz, vector<double>& v) { //ver si se pierde precisión con la división
  int n = matriz.dimensionFilas();
  int m = matriz.dimensionColumnas();
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      if(j == m - 1){
        v[i] = (v[i] + matriz[i][j]) / (float)m; //calculo promedio final
      } else {
        v[i] = (v[i] + matriz[i][j]); //acumulo
      } 
    }
  }
}

void calcular_matriz(Matriz& matriz) {
  int n = matriz.dimensionFilas();
  int m = matriz.dimensionColumnas();

  vector<double> media(n);
  calcular_media(matriz, media);
  double divisor = sqrt ((double)n - 1.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      matriz[i][j] = (matriz[i][j] - media[i]) / divisor;
    }
  }
  /* Cuando termina este paso tengo la matriz de covarianza en matriz */
  //ahora tengo que diagonalizar matriz
  Matriz *base_ortonormal = new Matriz(n, n);
  vector<double> aux(base_ortonormal->dimensionColumnas());  
  double autovalor;

  for (int i = 0; i < m; ++i) {
    Matriz::cargarVector(aux);
    autovalor = metodoPotencia(matriz, aux); //calculo el i-ésimo autovalor
    for (int k = 0; k < base_ortonormal->dimensionColumnas(); ++k){
      (*base_ortonormal)[i][k] = aux[k]; //completo matriz con base ortonormal
    }
    /* Deflación */
    Matriz auxiliar = Matriz::multiplicarVectoresDameMatriz(aux, aux);
    matriz- (auxiliar * autovalor);
  } 
  /* Tengo en matriz_ortonormal la matriz con base de autovectores de matriz */
  
}
