#include "clases/Matriz.h"
//#include "metodos/PLSDA.cpp"
//#include "metodos/PLSDATest.cpp"
#include "metodos/pca.cpp"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream> 
#include <algorithm>
#include <map>
#include <queue>
#include <sys/time.h>
#include <limits>

timeval sstart, eend;
double acum = 0;
double acum2 = 0;

void init_time() {
  gettimeofday(&sstart, NULL);
}

double get_time() {
  gettimeofday(&eend, NULL);
  return (1000000*(eend.tv_sec-sstart.tv_sec) + (eend.tv_usec-sstart.tv_usec))/1000000.0;
}

/* Auxiliares estadística */

void calcularEstadisticas(Matriz& imagenes) {
  std::vector<int> tp(10, 0);
  std::vector<int> fp(10, 0);
  std::vector<int> fn(10, 0);

  vector<double> precisiones(10, 0);
  vector<double> recalls(10, 0);

  int filas = imagenes.dimensionFilas();

  for (size_t i = 0; i < filas; ++i) {

    int estimacion = imagenes.dameEstimacion(i);
    int etiqueta = imagenes.dameEtiqueta(i);
    
    if (estimacion == etiqueta) { 
      tp[etiqueta] += 1;  // pertenecia a la clase "etiqueta" y lo estime bien. Es un verdadero positivo para la clase "etiqueta"
      
    } else { 
      fp[estimacion] += 1; // Supuse que era de la clase "estimacion" pero era de la clase "etiqueta". Es un falso positivo para la clase "estimacion"
      fn[etiqueta] += 1; // Pertenecia a la clase "etiqueta" pero me identificaron con la clase "estimacion". Es un falso negativo para la clase "etiqueta"
    }

  }

  //std::cout << "-----------------------------------------" << std::endl;
  //std::cout << "Cálculo precision: " << std::endl;
  for (int i = 0; i < 10; ++i) {
    double precision = (double)tp[i] / (double)(tp[i] + fp[i]);
    precisiones[i] = precision;

    double recall = (double)tp[i] / (double)(tp[i] + fn[i]);
    recalls[i] = recall;

    //std::cout << "Precisión clase " << i << ": " << precision << std::endl; 
  }

  std::cout << "-----------------------------------------" << std::endl;
  std::cout << "Cálculo precision: " << std::endl;


  

}

/* Auxiliares */

double dameNorma(std::vector<double>& x) {
  double norm = 0;
  int length = x.size();
  for (int i = 0; i < length; ++i)
    norm += x[i] * x[i];
    
  return sqrt(norm);
}

int indiceMinimo(std::vector<double>& x) {
  
  std::vector<double>::iterator itMin = std::min_element(x.begin(), x.end()); 
  int posMin = std::distance(x.begin(), itMin);

  x[posMin] = std::numeric_limits<double>::max(); // Tengo que hacer esto porque calculo k minimos. No se me ocurre otra manera ahora


  return posMin;
}

int indiceMaximo(std::vector<int>& x) {

  std::vector<int>::iterator itMax = std::max_element(x.begin(), x.end()); 
  int posMax = std::distance(x.begin(), itMax);

  return posMax;
}

/* KNN */
    
int dameEtiquetaEstimada(Matriz& imagenesTrain, std::vector<double>& imagen, int vecinos) {

  int filas = imagenesTrain.dimensionFilas();
  int columnas = imagenesTrain.dimensionColumnas();

  std::vector<double> x(columnas, 0);
  std::vector<double> y(filas, 0);

  for (int i = 0; i < filas; ++i) {

    for (int j = 0; j < columnas; ++j) 
      x[j] = imagen[j] - imagenesTrain[i][j];

    y[i] = dameNorma(x);
  }

  // El vector y tiene en la i-esima posicion la norma |imagen - y[i]|_2

  std::vector<int> labels(10, 0);

  while (vecinos > 0) {
    int i = indiceMinimo(y); // Me fijo cual es la imagen que minimiza la norma en cada iteracion
    int etiqueta = imagenesTrain.dameEtiqueta(i); // Me fijo cual es la etiqueta de dicho minimo
    labels[etiqueta] += 1;

    vecinos--;
  }

  return indiceMaximo(labels); // Devuelvo el "más votado"
 
}

int KNN(Matriz& imagenesTrain, Matriz& imagenesTest, int vecinos) {

  int filas = imagenesTest.dimensionFilas();
  int hitRate = 0;

  //std::cout << "----------------------------------------------------------" << std::endl;
  //imagenesTrain.mostrar2();
  //std::cout << "----------------------------------------------------------" << std::endl;
  //imagenesTest.mostrar2();
  //std::cout << "----------------------------------------------------------" << std::endl;
  
  for (int i = 0; i < filas; ++i) {

    int etiqueta = dameEtiquetaEstimada(imagenesTrain, imagenesTest[i], vecinos); // Le asigna a qué número pertenece la i-ésima imagen de imagenesTest
    imagenesTest.estimar(i, etiqueta); // En matriz.estimar tengo lo que supongo que es la imagen. En matriz.etiqueta tengo lo que de verdad es

 /*   
    std::cout << "la etiqueta de la imagen " << i << " es " << imagenesTest.dameEtiqueta(i) << std::endl;
    std::cout << "su estimación fue: " << imagenesTest.dameEstimacion(i) << std::endl;

    std::cout << "etiqueta: " << imagenesTest.dameEtiqueta(i) << std::endl;
    std::cout << "estimación: " << imagenesTest.dameEstimacion(i) << std::endl;
  */
    if (etiqueta == imagenesTest.dameEtiqueta(i))
      hitRate++;
  }

  std::cout << "cantidad de imagenes Test" << filas << std::endl;
  std::cout << "cantidad de aciertos (Hit Rate): " << hitRate << std::endl;

  calcularEstadisticas(imagenesTest);

  return hitRate; // Esto no sé si es necesario aún

}

int evaluarTests(std::string fileTestData, std::string fileTestResult, int method) {

  init_time();
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

  std::ifstream fileTest (test.c_str()); // Hasta aca sólo instancie variables

  while (getline (fileData, lineData)) { // Pido las K lineas

    std::ifstream fileTrain (train.c_str());

    std::istringstream issData(lineData);

    std::map<int, bool> isTrain;
    std::map<int, bool>::iterator itTrain = isTrain.begin();
  
    bool train;
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

    // Ahora tengo el tamaño de la matrix para poder instanciarla
    //

    int tamImagen = 784;

    Matriz imagenesTrain(cantImagenesTrain, tamImagen);
    Matriz imagenesTrainPCA(cantImagenesTrain, componentes);
    Matriz imagenesTrainPLSDA(cantImagenesTrain, dimensiones);

    Matriz imagenesTest(cantImagenesTest, tamImagen);
    Matriz imagenesTestPCA(cantImagenesTest, componentes);
    Matriz imagenesTestPLSDA(cantImagenesTest, dimensiones);

    int h = 0; // Marca el índice de imagenesTrain
    int r = 0; // Marca el índice de imagenesTest
    getline (fileTrain, lineTrain); // Descarto la primer línea que sólo tiene strings "label, pixel0, ..."

    for (int k = 0; k < cantImagenesTotales; ++k) {
    
      getline (fileTrain, lineTrain);

      std::istringstream issTrain(lineTrain);

      if (isTrain[k]) {

        int j = 0;
        int a = 0;
        std::string s;

        while (getline(issTrain, s, ',')) {
          if (j == 0) {
            imagenesTrain.etiquetar(h, atoi(s.c_str()));
            imagenesTrainPLSDA.etiquetar(h, atoi(s.c_str()));
            imagenesTrainPCA.etiquetar(h, atoi(s.c_str()));
            ++j;
          } else {
            imagenesTrain[h][a] = std::stod(s);
            ++a;
          }
        }

        ++h;

      } else {

        int j = 0;
        int a = 0;
        std::string s;

        while (getline(issTrain, s, ',')) {
          if (j == 0) {
            imagenesTest.etiquetar(r, atoi(s.c_str()));
            imagenesTestPLSDA.etiquetar(r, atoi(s.c_str()));
            imagenesTestPCA.etiquetar(r, atoi(s.c_str()));
            ++j;
          } else {
            imagenesTest[r][a] = std::stod(s);
            ++a;
          }
        }

        ++r;
      }
    }


    // imagenesTrain tiene en cada fila una imagen para la parte de train Y su correspondiente etiqueta (de 0 a 9)
    // imagenesTest tiene en cada fila una imagen para la parte de test Y su correspondiente etiqueta (de 0 a 9)
    // imagenesTrainPSLA tendrá la transformación característica de PSLA
    // imagenesTrainPCA tendrá la transformación característica de PCA
    // isTrain tiene para cada imagen (de 0 a cantImagenesTotales) si es parte del train o no

    Matriz autovectoresPCA(componentes, tamImagen);
    PCAMethod(imagenesTrain, imagenesTrainPCA, autovectoresPCA, componentes, fileWrite);

    Matriz autovectoresPLSDA(dimensiones, tamImagen);
    //PLSDAMethod(imagenesTrain, imagenesTrainPLSDA, autovectoresPLSDA, dimensiones, fileWrite); // Por ahora le hardcodeo el segundo parametro TODO: ver como cambiarlo
    

    switch(method) {
      case 0: { // Método KNN

          KNN(imagenesTrain, imagenesTest, vecinos);

          break;
              
        }

      case 1: { // Método KNN+PCA

        Matriz wTras(tamImagen, componentes);
        autovectoresPLSDA.trasponer(wTras);

        imagenesTest.multiplicarMatrices(wTras, imagenesTestPCA);

        KNN(imagenesTrainPCA, imagenesTestPCA, vecinos);
              
          break;
        }

      case 2: { // Método KNN+PLS-DA

        Matriz wTras(tamImagen, dimensiones);
        autovectoresPLSDA.trasponer(wTras);

        imagenesTest.multiplicarMatrices(wTras, imagenesTestPLSDA);

        KNN(imagenesTrainPLSDA, imagenesTestPLSDA, vecinos);
        break;

        }

    }

    ++z;
    std::cout << z << std::endl; // SOlo uso esto para ver cuantas iteraciones de lineas de archivo hizo

  }

  acum += get_time();
  std::cout << std::fixed << acum << std::endl;

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
