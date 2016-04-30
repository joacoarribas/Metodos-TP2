#include "clases/Matriz.h"
#include "clases/MatrizSimetrica.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream> 
#include <algorithm>
#include <map>
#include <queue>
#include <sys/time.h>

typedef std::priority_queue<int, std::vector<int>, std::greater<int> > min_heap; 

double dameNorma(std::vector<double> x) {
  double norm = 0;
  for (size_t i = 0; i < 784; ++i)
    norm += x[i] * x[i];
    
  return sqrt(norm);
}

int indiceMinimo(std::vector<double> x) {
  int min = 99999;
  int length = x.size();
  int j = 0;

  for (size_t i = 0; i < length; ++i) {
    if (x[i] < min) {
      min = x[i];
      j = i;
    } 
  }

  return j;

}

int indiceMaximo(std::vector<int> x) {
  int max = 0;
  int length = x.size();
  int j = 0;

  for (size_t i = 0; i < length; ++i) {
    if (x[i] < max) {
      max = x[i];
      j = i;
    } 
  }

  return j;

}
    
char dameEtiqueta(Matriz& imagenesTrain, std::vector<double>& imagen, int vecinos) {

  int filas = imagenesTrain.dimensionFilas();
  int columnas = imagenesTrain.dimensionColumnas();

  std::vector<double> x(columnas);
  std::vector<double> y(filas);

  for (size_t i = 0; i < filas; ++i) {

    for (size_t j = 0; j < columnas; ++j) 
      x[i] = imagen[j] - imagenesTrain[i][j];

    y[i] = dameNorma(x);
      // quiero sacar ccantidad de vecinos etiquetas (de 0 a 9). Sacar el maximo de ahí.
      
    }

  std::vector<int> labels(9, 0);

  while (vecinos > 0) {
    int i = indiceMinimo(y);
    char etiqueta = imagenesTrain.dameEtiqueta(i);
    labels[etiqueta]++;

    vecinos--;
  }

  return indiceMaximo(labels);
 
}

void KNN(Matriz& imagenesTrain, Matriz& imagenesTest, int vecinos) {

  int filas = imagenesTest.dimensionFilas();
  
  for (size_t i = 0; i < filas; ++i) {

    char etiqueta = dameEtiqueta(imagenesTrain, imagenesTest[i], vecinos); // Le asigna a qué número pertenece la i-ésima imagen de imagenesTest
    imagenesTest.etiquetar(i, etiqueta);
    
  }

}

int evaluarTests(std::string fileTestData, std::string fileTestResult, int method) {
  std::string lineData;
  std::string lineTest;
  std::string lineTrain;
  std::ifstream fileData (fileTestData.c_str());
  std::ofstream fileWrite (fileTestResult.c_str());
  int z = 0;

  getline(fileData, lineData); // Pido primer línea para instanciar todo
  std::istringstream issData(lineData);

  std::string path;
  std::string train;
  std::string test;
  int vecinos;
  int componentes;
  int dimensiones;
  int particiones;

  issData >> path;
  issData >> vecinos;
  issData >> componentes;
  issData >> dimensiones;
  issData >> particiones;

  train = path.append("train.csv"); 
  test = path.append("test.csv"); 

  std::ifstream fileTrain (train.c_str());
  std::ifstream fileTest (test.c_str());

  while (getline (fileData, lineData)) { // Pido las k lineas

    std::istringstream issData(lineData);

    std::map<int, bool> isTrain;
    std::map<int, bool>::iterator itTrain = isTrain.begin();
  
    char train;
    int cantImagenesTrain = 0;
    int cantImagenesTest = 0;
    int cantImagenesTotales = 0;

    int i = 0;
    while (issData >> train) { // Me fijo qué y cuántas imagenes van a ser utilizadas para el train. Idem test

      if (train) {
        cantImagenesTrain++;
        isTrain.insert(itTrain, std::pair<int, bool> (i, true));
      } else {
        cantImagenesTest++;
        isTrain.insert(itTrain, std::pair<int, bool> (i, false));
      }

      cantImagenesTotales++;

      ++itTrain;
      ++i;
    }

    // Ahora tengo el tamaño de la matrix  para poder instanciarla

    int tamImagen = 784;
    Matriz imagenesTrain(cantImagenesTrain, tamImagen);
    Matriz imagenesTest(cantImagenesTest, tamImagen);

    int h = 0; // Marca el índice de imagenesTrain
    int r = 0; // Marca el índice de imagenesTest
    getline (fileTrain, lineTrain); // Descarto la primer línea que sólo tiene strings

    for (int k = 0; k < cantImagenesTotales; ++k) {
    
      getline (fileTrain, lineTrain);

      std::istringstream issTrain(lineTrain);

      if (isTrain[k]) {

        char etiqueta;
        issTrain >> etiqueta;
        imagenesTrain.etiquetar(h, etiqueta);

        double pixel;
        int j = 0;

        while (issTrain >> pixel) {
          imagenesTrain[h][j] = pixel;
          ++j;
        }

        ++h;

      } else {

        char etiqueta;
        issTrain >> etiqueta;
        imagenesTrain.etiquetar(r, etiqueta);

        double pixel;
        int j = 0;

        while (issTrain >> pixel) {
          imagenesTrain[r][j] = pixel;
          ++j;
        }

        ++r;
      }
    }

    // imagenesTrain tiene en cada fila una imagen para la parte de train Y su correspondiente etiqueta (de 0 a 9)
    // imagenesTest tiene en cada fila una imagen para la parte de test Y su correspondiente etiqueta (de 0 a 9)
    // isTrain tiene para cada imagen (de 0 a cantImagenesTotales) si es parte del train o no
    
    switch(method) {
      case 0: { // Método KNN

          KNN(imagenesTrain, imagenesTest, vecinos);
              
        }

      case 1: { // Método KNN+PCA

              
        }

      case 2: { // Método KNN+PLS-DA

              
        }

    }


    ++z;
    std::cout << z << std::endl; // SOlo uso esto para ver cuantas iteraciones de lineas de archivo hizo

  }

  return 0;
}

int main(int argc, char** argv) {
  std::string fileTestData(argv[1]);
  std::string fileTestResult(argv[2]);
  int method(atoi(argv[3]));
  // Recibo por parametro tres archivos
  // El primero del cual leo los datos a evaluar
  // El segundo en el cual evaluo si los resultados fueron correctos
  // El tercero el método a realizar (0 KNN, 1 PCC+KNN, 2 PLS-DA+KNN)

  evaluarTests(fileTestData, fileTestResult, method);

  return 0;
}
