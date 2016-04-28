#include "clases/Matriz.h"
#include "clases/MatrizSimetrica.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <sys/time.h>

int evaluarTests(std::string fileTestData, std::string fileTestResult, int method) {
  std::string lineData;
  std::string lineTest;
  std::string lineTrain;
  std::ifstream fileData (fileTestData.c_str());
  std::ofstream fileWrite (fileTestResult.c_str());

  while (getlineData (fileData, lineData)) {

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

    int tamImagen = 784;
    int cantImagenes = 42000

    std::vector<std::list<int> > etiquetas(9);
    Matriz imagenes(cantImagen, tamImagen);

    for (int k = 0; k < cantImagenes; ++k) {
    
      getline (fileData, line);

      std::istringstream iss(line);

      std::string fecha;
      int equipo1;
      int equipo2;
      int golesEquipo1;
      int golesEquipo2;

      iss >> fecha;
      iss >> equipo1;
      iss >> golesEquipo1;
      iss >> equipo2;
      iss >> golesEquipo2;

      equipo1 -= 1;
      equipo2 -= 1;

      if (golesEquipo1 > golesEquipo2) { // Gano el equipo 1
        w[equipo1]++;
        l[equipo2]++;
      } else {
        w[equipo2]++;
        l[equipo1]++;
      }

      C[equipo1][equipo2]--;
      C[equipo2][equipo1]--;

    }

    std::vector<float> r(cantEquipos, 0);
    std::vector<float> y(cantEquipos, 0);

    for (int i = 0 ; i < cantEquipos ; ++i) {
      C[i][i] = 2.0 + w[i] + l[i];
      b[i] = 1.0 + (w[i] - l[i]) / 2.0; 
    }

    Matriz L(cantEquipos, cantEquipos);

    switch(method) {

    }

    for (int i = 0 ; i < cantEquipos ; ++i) {
      fileWrite << "equipo: " << i << " ranking: " << std::fixed << r[i] << std::endl; 
    }

    ++z;
    std::cout << z << std::endl;

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
