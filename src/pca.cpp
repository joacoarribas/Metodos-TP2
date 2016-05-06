#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "metodoPotencia.cpp"

vector<double>& calcular_media(Matriz& m) {
  double n = (double)m.dimensionFilas();
  vector<double> aux(m.dimensionFilas());
  for (int i = 0; i < m.dimensionFilas(); ++i) {
    for (int j = 0; j < m.dimensionColumnas(); ++j) {
      if(j == m.dimensionColumnas() - 1){
        aux[i] = (aux[i] + m[i][j]) / m.dimensionColumnas(); //no se si se pierde precisión así
      } else {
        aux[i] = (aux[i] + m[i][j]); //no se si se pierde precisión así
      } 
    }
  }
  return aux;
}

void calcular_matriz(Matriz& m) {
  vector<double> media = calcular_media(m);
  double n = sqrt ((double)m.dimensionFilas() - 1.0);
  for (int i = 0; i < m.dimensionFilas(); ++i) {
    for (int j = 0; j < m.dimensionColumnas(); ++j) {
      m[i][j] = (m[i][j] - media[i]) / n;
    }
  }
  /* Cuando termina este paso tengo la matriz de covarianza */
  //ahora tengo que diagonalizar m
  Matriz *base_ortonormal = new Matriz(m.dimensionFilas(), m.dimensionFilas());
  vector<double> aux(base_ortonormal->dimensionColumnas());  
  double autovalor;
  for (int i = 0; i < m.dimensionColumnas(); ++i) {
    Matriz::cargarVector(aux);
    autovalor = metodoPotencia(m, aux);
    for (int k = 0; k < base_ortonormal->dimensionColumnas(); ++k){
      (*base_ortonormal)[i][k] = aux[k];
    }
    Matriz auxiliar = Matriz::multiplicarVectoresDameMatriz(aux, aux);
    m - (auxiliar * autovalor);
  } 
  
}
