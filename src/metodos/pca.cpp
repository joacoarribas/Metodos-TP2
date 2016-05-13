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
}

void calcular_matriz_covarianza(Matriz& matriz, Matriz& res){
  int n = matriz.dimensionFilas();
  int m = matriz.dimensionColumnas();

  vector<double> media(m, 0);
  
  calcular_media(matriz, media);
  
  double divisor = sqrt ((double)n - 1.0);

  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i) {
      res[i][j] = (matriz[i][j] - media[j]) / divisor; //escribo resultado en la matriz
    }
  }
}

void calcular_base_ortonormal(Matriz& matriz, Matriz& matriz_ortonormal, int alfa){ //deja en matriz_ortonormal una matriz de alfa columnas
  double autovalor;
  int m = matriz_ortonormal.dimensionColumnas();
  vector<double> aux(m);

  for (int i = 0; i < alfa; ++i) { //repito alfa veces (hay que experimentar con dicho valor)
    Matriz::cargarVector(aux);

    autovalor = metodoPotencia(matriz, aux); //calculo el i-ésimo autovalor, en aux queda el autovector
    for (int k = 0; k < m; ++k){
      matriz_ortonormal[i][k] = aux[k]; //completo matriz con base ortonormal
    }
    /* Deflación */
    Matriz auxiliar(m, m);
    auxiliar.multiplicarVectoresDameMatriz(aux, aux);
    auxiliar.multiplicarEscalar(autovalor);
    matriz.menos(auxiliar);
  } 
  /* Tengo en matriz_ortonormal la matriz con base de autovectores de matriz */
}

void PCA(Matriz& matriz, Matriz& res, int alfa){
  int n = matriz.dimensionFilas();
  int m = matriz.dimensionColumnas();

  Matriz X(n, m);
  std::cout << "1" << std::endl;
  calcular_matriz_covarianza(matriz, X);
  std::cout << "2" << std::endl;

  Matriz Xt(m, n);
  X.trasponer(Xt);
  
  Matriz m_covarianza(m, m); //revisar dimensiones

  Xt.multiplicarMatrices(X, m_covarianza); //crear matriz covarianza

  Matriz m_ortonormal(alfa, m);

  /* Reducción de la dimensión */
  std::cout << "3" << std::endl;
  calcular_base_ortonormal(m_covarianza, m_ortonormal, alfa); //deja en matriz_ortonormal una matriz de alfa filas
  std::cout << "4" << std::endl;
  /* Transformación característica */
  Matriz m_ortonormal_traspuesta(m, alfa);
  std::cout << "5" << std::endl;
  m_ortonormal.trasponer(m_ortonormal_traspuesta);
  std::cout << "6" << std::endl;

  /* Transformación característica */
  matriz.multiplicarMatrices(m_ortonormal_traspuesta, res); // Chequear esto en caso de que no ande
  std::cout << "7" << std::endl;

}
