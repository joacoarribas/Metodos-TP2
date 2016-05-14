#include "Matriz.h"

Matriz::Matriz(int filas, int columnas) {
  this->filas    = filas;
  this->columnas = columnas;

  vector< vector<double> > m(filas, vector<double>(columnas));
  this->matriz = m;

  vector<int> v(filas, 0);
  this->etiquetas = v;
  this->estimaciones = v;

}

Matriz::~Matriz() {}

int Matriz::dameEtiqueta(int i) {
  return this->etiquetas[i];
}

int Matriz::dameEstimacion(int i) {
  return this->estimaciones[i];
}

int Matriz::dimensionFilas() {
  return this->filas;
}

int Matriz::dimensionColumnas(){
  return this->columnas;
}

void Matriz::etiquetar(int i, int etiqueta){
  this->etiquetas[i] = etiqueta;
}

void Matriz::estimar(int i, int etiqueta){
  this->estimaciones[i] = etiqueta;
}

void Matriz::trasponer(Matriz& traspuesta) {
  assert(this->filas == traspuesta.dimensionColumnas() && this->columnas == traspuesta.dimensionFilas());

  for (int f = 0; f < this->filas; ++f)
    for (int c = 0; c < this->columnas; ++c)
      traspuesta[c][f] = this->matriz[f][c];

}

	/*******************************
	 *          Operadores         *
	 *******************************/ 

void Matriz::multiplicarEscalar(double i) {

  for (int f = 0; f < this->filas; ++f)
    for (int c = 0; c < this->columnas; ++c)
      this->matriz[f][c] *= i;
  
}

void Matriz::multiplicarVectorDer(vector<double>& x, vector<double>& y) {
  assert(this->columnas == x.size());
  assert(this->filas == y.size());

  for (int f = 0; f < this->filas; ++f)
    for (int k = 0; k < this->columnas; ++k)
      y[f] += this->matriz[f][k] * x[k];

}

void Matriz::multiplicarVectorIzq(vector<double>& x, vector<double>& y) {
  assert(this->columnas == y.size());
  assert(this->filas == x.size());

  for (int k = 0; k < this->columnas; ++k)
    for (int f = 0; f < this->filas; ++f)
      y[k] += this->matriz[f][k] * x[f];

}

void Matriz::multiplicarMatrices(Matriz& a, Matriz& b) {
  assert(this->columnas == a.dimensionFilas());
  assert(b.dimensionFilas() == this->filas && b.dimensionColumnas() == a.dimensionColumnas());

  int col = a.dimensionColumnas();
  for (int f = 0; f < this->filas; ++f)
    for (int c = 0; c < col; ++c)
      for (int k = 0; k < this->columnas; ++k)
        b[f][c] += this->matriz[f][k] * a.matriz[k][c];
      
}

void Matriz::menos(Matriz& m) {
  assert (m.dimensionFilas() == this->filas && m.dimensionColumnas() == this->columnas);
  
    for (int f = 0; f < this->filas; ++f)
      for (int c = 0; c < this->columnas; ++c)
        this->matriz[f][c] = this->matriz[f][c] - m.matriz[f][c];
  
}

void Matriz::mas(Matriz& m) {
  assert (m.dimensionFilas() == this->filas && m.dimensionColumnas() == this->columnas);
  
    for (int f = 0; f < this->filas; ++f)
      for (int c = 0; c < this->columnas; ++c)
        this->matriz[f][c] = this->matriz[f][c] + m.matriz[f][c];
  
}

vector<double>& Matriz::operator [] (int fila) {
      return this->matriz[fila];
}

void Matriz::randomizar(int semilla) {
  srand( semilla );

  for (int f = 0; f < filas; ++f) {
    for (int c = 0; c < columnas; ++c) {
      matriz[f][c] = rand() % 1000 + 1;
    }	
  }	
}

void Matriz::mostrar() {
  for (int f = 0; f < this->filas; ++f) {
    cout << "  ";
    for (int c = 0; c < 6; ++c) {
      if (matriz[f][c] >= 100) {
        cout << "| " << matriz[f][c] << " |" << " ";
      } else {
        if (matriz[f][c] >= 10) {
          cout << "|  " << matriz[f][c] << " |" << " ";
        } else {
          cout << "|  " << matriz[f][c] << "  |" << " ";
        }
      }
    }

  cout << endl;
  }
}

void Matriz::mostrar2() {
  for (int f = 4; f < 20; ++f) {
    cout << "  ";
    for (int c = 0; c < 7; ++c) {
      if (matriz[f][c] >= 100) {
        cout << "| " << matriz[f][c] << " |" << " ";
      } else {
        if (matriz[f][c] >= 10) {
          cout << "|  " << matriz[f][c] << " |" << " ";
        } else {
          cout << "|  " << matriz[f][c] << "  |" << " ";
        }
      }
    }

  cout << endl;
  cout << endl;
  cout << endl;
  }
}

void Matriz::cargarVector(std::vector<double> &x) {
  //para cargar el vector
  srand (time(NULL));

  for (int i=0; i<x.size(); ++i) {
    x[i] = (rand() % 10) + 1;
  }
}

void Matriz::cerearVector(std::vector<double> &x) {
  int n = x.size();

  for (int i = 0; i < n; ++i)
    x[i] = 0;

}

void Matriz::multiplicarVectoresDameMatriz(std::vector<double>& a, std::vector<double>& b) {
  assert(this->filas == a.size() && this->columnas == b.size());
  int columnas = b.size();
  int filas = a.size();

  for (int i=0; i < filas ; ++i)
    for (int j=0; j < columnas ; ++j)
      this->matriz[i][j] = a[i] * b[j];

}

double Matriz::multiplicarVectoresDameValor(std::vector<double>& a, std::vector<double>& b) {
  assert(a.size() == b.size());

  double result = 0;
  int length = a.size();

  for (int i = 0; i < length; ++i) {
    result += a[i] * b[i];
  }

  return result;
}

	/*******************************
	 *     Funciones privadas      *
	 *******************************/ 

void Matriz::clean() {
  for (int f = 0; f < filas; ++f) {
    for (int c = 0; c < columnas; ++c) {
      matriz[f][c] = 0;
    }
  }
}
  /*
  std::cout << "Matriz this:" << std::endl;
  std::cout << "filas: " << this->filas << std::endl;
  std::cout << "columnas: " << this->columnas << std::endl;

  std::cout << "_------------------------" << std::endl;

  std::cout << "Matriz parametro:" << std::endl;
  std::cout << "filas: " << a.dimensionFilas() << std::endl;
  std::cout << "columnas: " << a.dimensionColumnas() << std::endl;

  std::cout << "_------------------------" << std::endl;

  std::cout << "Matriz nueva:" << std::endl;
  std::cout << "filas: " << b.dimensionFilas() << std::endl;
  std::cout << "columnas: " << b.dimensionColumnas() << std::endl;
*/

