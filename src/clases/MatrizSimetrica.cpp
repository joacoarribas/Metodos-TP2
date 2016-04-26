#include "MatrizSimetrica.h"

MatrizSimetrica::MatrizSimetrica(int filas, int columnas) {
  this->filas    = filas;
  this->columnas = columnas;
  vector< vector<float> > m(filas, vector<float>(filas));
  this->matriz = m;
  this->transpuesta = false;
  this->triangular  = false;

  for (int i = 0; i < filas; ++i) {
    vector<float> fila(i+1); // <- para la matriz simetrica
    matriz[i] = fila;
  }
}

int MatrizSimetrica::dimensionFilas() {
  return this->filas;
}

int MatrizSimetrica::dimensionColumnas() {
  return this->columnas;
}

    /*******************************
     *          Operadores         *
     *******************************/ 

float MatrizSimetrica::get(int f, int c) {
  if (transpuesta) {
    int aux = c;
    c = f;
    f = aux;
  } 


  if (c <= f) {
    return this->matriz[f][c];
  } else {
    if (this->esTriangular()) {
      return 0;
    } else {
      return this->matriz[c][f];
    }
  }
}

void MatrizSimetrica::set(int f, int c, float valor) {
  if (transpuesta) {
    int aux = c;
    c = f;
    f = c;
  } 
  
  if (c <= f) {
    matriz [f][c] = valor;
  }
}

MatrizSimetrica& MatrizSimetrica::operator * (int i) {
  for (int f = 0; f < this->filas; ++f) {
    for (int c = 0; c < this->columnas; ++c) {
      this->set(f, c, this->get(f,c)*i);
    }
  }
  
  return *this;
}

Matriz& MatrizSimetrica::operator * (Matriz& m) {
  assert(this->columnas == m.dimensionFilas());

  Matriz * nueva = new Matriz(this->columnas, m.dimensionFilas());

  for (int f = 0; f < this->filas; ++f) {
    for (int c = 0; c < this->columnas; ++c) {
      for (int k = 0; k < this->columnas; ++k) {
        (*nueva)[f][c] += this->get(f,k) + m[k][c];
      }
    }
  }

  return *nueva;
}

Matriz& MatrizSimetrica::operator * (MatrizSimetrica& m) {
  assert(this->columnas == m.dimensionFilas());

  Matriz * nueva = new Matriz(this->columnas, m.dimensionFilas());

  for (int f = 0; f < this->filas; ++f) {
    for (int c = 0; c < this->columnas; ++c) {
      for (int k = 0; k < this->columnas; ++k) {
        (*nueva)[f][c] += this->get(f,k) * m.get(k,c);
      }
    }
  }

  return *nueva;
}

void MatrizSimetrica::randomizar(int semilla) {
  srand( semilla );

  for (int f = 0; f < filas; ++f) {
    for (int c = 0; c <= f; ++c) {
      this->set(f,c, rand() % 1000 + 1);
    }   
  } 
}

void MatrizSimetrica::mostrar() {
  for (int f = 0; f < filas; ++f) {
    cout << "  ";
    //for (int c = 0; c < this->matriz[f].size(); ++c) {
    for (int c = 0; c < columnas; ++c) {
      if (this->get(f,c) < 0) {
        cout << "| " << this->get(f,c) << "  |" << " ";
      } else {
        cout << "|  " << this->get(f,c) << "  |" << " ";
      }
    }

  cout << endl;
  }
}

void MatrizSimetrica::setTriangular(bool esTriangular) {
  triangular = esTriangular;
}

bool MatrizSimetrica::esTriangular() {
  return this->triangular;
}

void MatrizSimetrica::transponer() {
  transpuesta = !transpuesta;
}

    /*******************************
     *     Funciones privadas      *
     *******************************/ 

void MatrizSimetrica::clean() {
  for (int f = 0; f < filas; ++f) {
    for (int c = 0; c <= f; ++c) {
      matriz[f][c] = 0;
    }
  }
}
