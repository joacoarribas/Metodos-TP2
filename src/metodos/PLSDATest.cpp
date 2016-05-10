#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "metodoPotencia.cpp"

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

void PLSDAMethod(Matriz& imagenes, Matriz& imagenesTransformadas, int dimensiones) {
	//revizar imagenes de 1...n, cada uno es vec de m
	long n = imagenes.dimensionFilas();
	long m = imagenes.dimensionColumnas();
	double raizN = sqrt((double)n-1.0);

	Matriz X(n, m);
	Matriz Y(n, 10);

	//calcular mu = promedio de la suma de todos los vectores
  std::vector<double> promedios(n, 0);
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i) {
			promedios[j] += imagenes[i][j];
		}

		promedios[j] /= (double)n; 
	}

	//X matriz de nXm donde cada fila es traspuesta(x_i − μ)/raiz(n − 1)
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i) {
			X[i][j] = (imagenes[i][j] - promedios[j]) / raizN;
		}
	}

	//defino preY de nX10 que tiene 1 en preY_i,j si la i-esima muestra de la base corresponde al digito j-1
  //-1 en caso contrario
	Matriz preY(n, 10);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < 10; ++j) {
			if (imagenes.dameEtiqueta(i) == j) {
				preY[i][j] = 1;
			} else {
				preY[i][j] = -1;
			}
		}
  }

	//defino promY de nX1 con promY_i = promedio de la fila i-esima de preY
  std::vector<double> promY(n, 0);
  for (int j = 0; j < 10; ++j) {
    for (int i = 0; i < n; ++i) {
			promY[j] += preY[i][j];
		}
		promY[j] /= (double)n; //esto puede traer error
	}

	//defino Y como preY_i − promY_i / raiz(n-1)
  for (int j = 0; j < 10; ++j) {
    for (int i = 0; i < n; ++i) {
			Y[i][j] = (preY[i][j] - promY[j]) / raizN;
		}
	}

  //Matriz w(m, dimensiones); // Solo necesito dimensiones cantidad de autovectores
  Matriz w(dimensiones, m); // Solo necesito dimensiones cantidad de autovectores

  // Tengo dimensiones cantidad de autovectores de tamaño m

	for(int i=0; i<dimensiones; ++i) {

		Matriz Xtras(m, n);
    X.trasponer(Xtras);

		Matriz Ytras(10,n);
    Y.trasponer(Ytras);

		Matriz Mi(m, m);
		
    if (true) { // Hago esto para Xtranss e Ytranss mueran despues del scope, no las necesito más
      Matriz Xtrass(m, 10);
      Matriz Ytrass(10, m);

      Xtras.multiplicarMatrices(Y, Xtrass); // Xtrass = X^t * Y
      Ytras.multiplicarMatrices(X, Ytrass); // Ytrass = Y^t * X

      Xtrass.multiplicarMatrices(Ytrass, Mi); // Xtrass * Ytrass = X^t * Y * Y^t * X
    }

		//calcular wi el autovector asociado al mayor autovalo de Mi
		std::vector<double>& wi = w[i]; // Si w = nxm entonces wi = m
		Matriz::cargarVector(wi);
    metodoPotencia(Mi, wi); //descarto el autovalor que vino, solo necesito el auvector en wi

		//normalizar wi con norma 2
		normalizar(wi);
		
		//ti = X * wi
		std::vector<double> ti(n, 0);
    X.multiplicarVectorDer(wi, ti); // x = nxm y wi = mx1 entonces ti = nx1; 
		
		//normalizar ti con norma 2
		normalizar(ti);
		
		//	actualizar X como X − ti^t * ti * X;
    std::vector<double> tiX(m, 0);
    X.multiplicarVectorIzq(ti, tiX); // ti = 1xn y X = nxm entonces tiX = 1xm ; tiX = ti^t * X

    Matriz nX(n, m);
    nX.multiplicarVectoresDameMatriz(ti, tiX); // nX = ti * tiX = ti * ti^t * X

    X.menos(nX); // Termino de actualizar X = X - ti * ti^t * X

    std::vector<double> tiY(10, 0);
    Y.multiplicarVectorIzq(ti, tiY); // ti = 1xn y Y = nx10 entonces tiY = 1x10 ; tiY = ti^t * Y

    Matriz nY(n, 10);
    nY.multiplicarVectoresDameMatriz(ti, tiY); // nX = ti * tii = ti * ti^t * X

    Y.menos(nY); // Termino de actualizar X
    std::cout << "iteracion: " << i+1 << " de " << dimensiones << std::endl;
	}

  // Los autovectores estan como FILAS de la matriz w. Para multiplicarla por imagenes tengo que trasponerla
  
  Matriz wTras(m, dimensiones);

  w.trasponer(wTras);

  imagenes.multiplicarMatrices(wTras, imagenesTransformadas); // Esto es la transformación caracteritisca
}


/*
int main() {
	PLSDA * prueba = new PLSDA();

	return 0;
}*/
