#include "PLSDA.h"

PLSDA::PLSDA() {}

Matriz& PLSDA::PLSDAMethod(Matriz valores, int dimensiones) {
	//revizar valores de 1...n, cada uno es vec de m
	long n = valores.dimensionFilas();
	long m = valores.dimensionColumnas();
	double raizN = sqrt(n-1);
	Matriz X(n, m);
	Matriz Y(n, m);

	//calcular mu = promedio de la suma de todos los vectores
	double promedios[n];
	for (int i=0; i<n; ++i) {
		for (int j=0; j<m; ++j) {
			promedios[i] += valores[i][j];
		}
		promedios[i] /= m; //esto puede traer error
	}

	//X matriz de nXm donde cada fila es traspuesta(x_i − μ)/raiz(n − 1)
	for (int i=0; i<n; ++i) {
		for (int j=0; j<m; ++j) {
			X[i][j] = (valores[i][j] - promedios[i]) / raizN;
		}
	}

	//defino preY de nX10 que tiene 1 en preY_i,j si la i-esima muestra de la base corresponde al digito j-1
								 //-1 en caso contrario
	Matriz preY(n, 10);
	for (int i=0; i<n; ++i) {
		for (int j=0; j<10; ++j) {
			if (valores.dameEtiqueta(i) == j) {
				preY[i][j] = 1;
			} else {
				preY[i][j] = -1;
			}
		}
	}		

	//defino promY de nX1 con promY_i = promedio de la fila i-esima de preY
	double promY[n];
	for (int i=0; i<n; ++i) {
		for (int j=0; j<10; ++j) {
			promY[i] += valores[i][j];
		}
		promY[i] /= 10; //esto puede traer error
	}

	//defino Y como preY_i − promY_i / raiz(n-1)
	for (int i=0; i<n; ++i) {
		for (int j=0; j<m; ++j) {
			Y[i][j] = (preY[i][j] - promY[i]) / raizN;
		}
	}

	//genero las transformaciones w_i como :
	std::vector< std::vector<double> > w(n, std::vector<double>(n));
	for (int i = 0; i < n; ++i) {
	    std::vector<double> wi(n);
	    w[i] = wi;
	}

	for(int i=0; i<dimensiones; ++i) {
		Matriz Xtrans = X.transponer();
		Matriz Ytrans = Y.transponer();
		Matriz Mi     = (Xtrans * Y) * (Ytrans * X);

		//calcular wi el autovector asociado al mayor autovalo de Mi
		std::vector<double> &wi = w[i];
		cargarVector(wi);
		metodoPotencia(Mi, wi); //descarto el autovalor que vino, solo necesito el auvector en wi

		//normalizar wi con norma 2
		normalizar(wi);
		
		//ti = X * wi
		std::vector<double> ti = X * wi;
		
		//normalizar ti con norma 2
		normalizar(ti);
		
		//	actualizar X como X − ti^t * ti * X;
		X = X - (X * multiplicarVectores(ti, ti));

		//	actualizar Y como Y − ti^t * ti * Y;
		Y = Y - (Y * multiplicarVectores(ti, ti));
	}

	//aplico transformacion caracteristica (w) a las imagenes
	return transformacionCaracteristica(valores, w, n, dimensiones); 
}

double PLSDA::multiplicarVectores(std::vector<double> a, std::vector<double> b) {
	double result = 0;

	for (int i=0; i<a.size(); ++i) {
		result += a[i] * b[i];
	}

	return result;
}

void PLSDA::normalizar(std::vector<double> &x) {
	double norma2 = 0;
	for (int i=0; i<x.size(); ++i) {
		norma2 = norma2 + (x[i] * x[i]);
	}

	norma2 = sqrt(norma2);

	for (int i=0; i<x.size(); ++i) {
		x[i] = x[i] / norma2;
	}		
}

void PLSDA::cargarVector(std::vector<double> &x) {
	//para cargar el vector
	srand (time(NULL));

	for (int i=0; i<x.size(); ++i) {
		x[i] = rand() % 10;
	}
}

Matriz& PLSDA::transformacionCaracteristica(Matriz m, std::vector< std::vector<double> > w, int n, int dimensiones) {
	Matriz * resultado = new Matriz(n, dimensiones); // la nueva matriz usa el parametro de dimensiones

	for (int i=0; i<n; ++i) {
		for (int j=0; j<dimensiones; ++j) {
			for (int k=0; k<dimensiones; ++k){
				(*resultado)[i][j] = multiplicarVectores(w[i], m[k]); //multiplico cada imagen por la transformacion correspondiente
			}
		}
	}

	return *resultado;
}

/*
int main() {
	PLSDA * prueba = new PLSDA();

	return 0;
}*/