#include "Matriz.h"

Matriz::Matriz(int filas, int columnas) {
  this->filas    = filas;
  this->columnas = columnas;
  vector< vector<double> > m(filas, vector<double>(columnas));
  this->matriz = m;
  vector<char> v(filas, -1);
  this->etiquetas = v;
  this->estimaciones = v;

  for (int i = 0; i < filas; ++i) {
    vector<double> columna(filas);
    matriz[i] = columna;
  }

}

char Matriz::dameEtiqueta(int i) {
  return this->etiquetas[i];
}

char Matriz::dameEstimacion(int i) {
  return this->estimaciones[i];
}

int Matriz::dimensionFilas() {
  return this->filas;
}

int Matriz::dimensionColumnas(){
  return this->columnas;
}

void Matriz::etiquetar(int i, char etiqueta){
  this->etiquetas[i] = etiqueta;
}

void Matriz::estimar(int i, char etiqueta){
  this->estimaciones[i] = etiqueta;
}

void Matriz::transponer() {
  for (int f = 0; f < this->filas; ++f) {
    for (int c = f+1; c < columnas; ++c) {
      float aux = this->matriz[f][c];
      
      this->matriz[f][c] = this->matriz[c][f];
      this->matriz[c][f] = aux;
    }
  }
}

	/*******************************
	 *          Operadores         *
	 *******************************/ 

Matriz& Matriz::operator * (int i) {
  for (int f = 0; f < this->filas; ++f) {
    for (int c = 0; c < this->columnas; ++c) {
      matriz [f][c] *= i;
    }
  }
  
  return *this;
}

vector<double>& Matriz::operator * (vector<double>& v) {
  assert(this->columnas == v.size());

  vector<double> * nuevo = new vector<double>(v.size());

  for (int f = 0; f < this->filas; ++f) {
    for (int k = 0; k < this->columnas; ++k) {
      (*nuevo)[f] += this->matriz[f][k] * v[k];
    }
  }

  return *nuevo;
}

Matriz& Matriz::operator * (Matriz& m) {
  assert(this->columnas == m.dimensionFilas());

  Matriz * nueva = new Matriz(this->columnas, m.dimensionFilas());

  for (int f = 0; f < this->filas; ++f) {
    for (int c = 0; c < this->columnas; ++c) {
      for (int k = 0; k < this->columnas; ++k) {
        nueva->matriz[f][c] += this->matriz[f][k] * m.matriz[k][c];
      }
    }
  }

  return *nueva;
}

Matriz& Matriz::operator * (MatrizSimetrica& m) {
  assert(this->columnas == m.dimensionFilas());

  Matriz * nueva = new Matriz(this->columnas, m.dimensionFilas());

  for (int f = 0; f < this->filas; ++f) {
    for (int c = 0; c < this->columnas; ++c) {
      for (int k = 0; k < this->columnas; ++k) {
        nueva->matriz[f][c] += this->matriz[f][k] * m.get(k,c);
      }
    }
  }

  return *nueva;
}

Matriz& Matriz::operator + (Matriz& m) {
  if (m.dimensionFilas() == this->filas && m.dimensionColumnas() == this->columnas) {
  
    for (int f = 0; f < this->filas; ++f) {
      for (int c = 0; c < this->columnas; ++c) {
        this->matriz[f][c] = this->matriz[f][c] + m.matriz[f][c];
      }
    }
  
  }
  
  return *this;
}

Matriz& Matriz::operator - (Matriz& m) {
  if (m.dimensionFilas() == this->filas && m.dimensionColumnas() == this->columnas) {
  
    for (int f = 0; f < this->filas; ++f) {
      for (int c = 0; c < this->columnas; ++c) {
        this->matriz[f][c] = this->matriz[f][c] - m.matriz[f][c];
      }
    }
  
  }
  
  return *this;
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
  for (int f = 0; f < filas; ++f) {
    cout << "  ";
    for (int c = 0; c < columnas; ++c) {
      if (matriz[f][c] < 0) {
        cout << "| " << matriz[f][c] << "  |" << " ";
      } else {
        cout << "|  " << matriz[f][c] << "  |" << " ";
      }
    }

  cout << endl;
  }
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
