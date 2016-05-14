#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "metodoPotencia.cpp"

/* Funciones PCA*/

void normalizar(std::vector<double>& x) {
	double norma2 = 0;
  long length = x.size();

	for (int i=0; i<length; ++i) {
		norma2 = norma2 + (x[i] * x[i]);
	}

	norma2 = sqrt(norma2);

	for (int i=0; i<length; ++i) {
		x[i] = x[i] / norma2;
	}		
}

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

void calcular_matriz_X(Matriz& matriz, Matriz& res){
  int n = matriz.dimensionFilas();
  int m = matriz.dimensionColumnas();

  vector<double> media(m, 0.0);
  
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

  for (int i = 0; i < alfa; ++i) { //repito alfa veces (hay que experimentar con dicho valor)

    vector<double>& aux = matriz_ortonormal[i];

    Matriz::cargarVector(aux);

    matriz.mostrar2();

    autovalor = metodoPotencia(matriz, aux); //calculo el i-ésimo autovalor, en aux queda el autovector
    std::cout << "autovalor: " << i << " vale: " << std::scientific << autovalor << std::endl;

    /* Deflación */
    Matriz auxiliar(m, m);
    std::cout << "-----------------------------------" << std::endl;

    normalizar(aux);

    auxiliar.multiplicarVectoresDameMatriz(aux, aux);
    auxiliar.multiplicarEscalar(autovalor);
    matriz.menos(auxiliar);

    matriz.mostrar2();

    while (true) {}
  } 
  /* Tengo en matriz_ortonormal la matriz con base de autovectores de matriz */
}

void PCAMethod(Matriz& matriz, Matriz& res, Matriz& m_ortonormal, int alfa, std::ofstream& filewrite){
  int n = matriz.dimensionFilas();
  int m = matriz.dimensionColumnas();

  Matriz X(n, m);
  calcular_matriz_X(matriz, X);

  Matriz Xt(m, n);
  X.trasponer(Xt);

  Matriz m_covarianza(m, m); //revisar dimensiones

  std::cout << "antes de entrar a multiplicar" << std::endl;
  Xt.multiplicarMatrices(X, m_covarianza); //crear matriz covarianza
  std::cout << "after" << std::endl;

  /* Reducción de la dimensión */
  calcular_base_ortonormal(m_covarianza, m_ortonormal, alfa); //deja en matriz_ortonormal una matriz de alfa filas

  /* Transformación característica */
  Matriz m_ortonormal_traspuesta(m, alfa);
  m_ortonormal.trasponer(m_ortonormal_traspuesta);

  /* Transformación característica */
  matriz.multiplicarMatrices(m_ortonormal_traspuesta, res); // Chequear esto en caso de que no ande
}
