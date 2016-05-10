#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "metodos/metodoPotencia.cpp"

void calcular_media(Matriz& matriz, vector<double>& v) { //ver si se pierde precisión con la división
  int n = matriz.dimensionFilas();
  int m = matriz.dimensionColumnas();
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i) {
      v[j] = (v[j] + matriz[i][j]); //acumulo
    }
    v[j] = v[j] / (double)m; //calculo promedio final
  }
}

Matriz& calcular_matriz_covarianza(Matriz& matriz){
  int n = matriz.dimensionFilas();
  int m = matriz.dimensionColumnas();
  Matriz *res = new Matriz(n, n);

  vector<double> media(n);
  calcular_media(matriz, media);
  double divisor = sqrt ((double)n - 1.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      (*res)[i][j] = (matriz[i][j] - media[i]) / divisor; //escribo resultado en la matriz
    }
  }
  return *res;
}

void calcular_base_ortonormal(Matriz& m, Matriz& matriz_ortonormal, int alfa){ //deja en matriz_ortonormal una matriz de alfa columnas
  double autovalor;
  vector<double> aux(matriz_ortonormal.dimensionColumnas());  

  for (int i = 0; i < alfa; ++i) { //repito alfa veces (hay que experimentar con dicho valor)
    Matriz::cargarVector(aux);
    autovalor = metodoPotencia(m, aux); //calculo el i-ésimo autovalor, en aux queda el autovector
    for (int k = 0; k < matriz_ortonormal.dimensionColumnas(); ++k){
      matriz_ortonormal[i][k] = aux[k]; //completo matriz con base ortonormal
    }
    /* Deflación */
    Matriz auxiliar = Matriz::multiplicarVectoresDameMatriz(aux, aux);
    m - (auxiliar * autovalor);
  } 
  /* Tengo en matriz_ortonormal la matriz con base de autovectores de matriz */
}

Matriz& PCA(Matriz& matriz){
  int n = matriz.dimensionFilas();
  int m = matriz.dimensionColumnas();

  Matriz matriz_covarianza = calcular_matriz_covarianza(matriz);

  int alfa = 10; //alfa valor arbitrario
  Matriz *matriz_ortonormal = new Matriz(n, alfa); //revisar dimensiones
  /* Reducción de la dimensión */
  calcular_base_ortonormal(matriz, *matriz_ortonormal, alfa); //deja en matriz_ortonormal una matriz de alfa columnas
  /* Transformación característica */
  Matriz matriz_ortonormal_traspuesta = matriz_ortonormal->transponer();

  Matriz *res = new Matriz(alfa, n); //revisar dimensiones
  *res = matriz_ortonormal_traspuesta * matriz;
  return *res;
}

int main(){
  return 0;
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

