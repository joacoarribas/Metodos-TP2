#include "clases/Matriz.h"
#include "metodos/PLSDA.cpp"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream> 
#include <algorithm>
#include <map>
#include <queue>
#include <sys/time.h>

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

  x[posMin] = 99999999; // Tengo que hacer esto porque calculo k minimos. No se me ocurre otra manera ahora
  
  return posMin;
}

int indiceMaximo(std::vector<int>& x) {

  std::vector<int>::iterator itMax = std::max_element(x.begin(), x.end()); 
  int posMax = std::distance(x.begin(), itMax);

  return posMax;
}
    
int dameEtiqueta(Matriz& imagenesTrain, std::vector<double>& imagen, int vecinos) {

  int filas = imagenesTrain.dimensionFilas();
  int columnas = imagenesTrain.dimensionColumnas();

  std::vector<double> x(columnas);
  std::vector<double> y(filas);

  for (int i = 0; i < filas; ++i) {

    for (int j = 0; j < columnas; ++j) 
      x[j] = imagen[j] - imagenesTrain[i][j];

    y[i] = dameNorma(x);
//    std::cout << "norma: " << y[i] << std::endl;
      // quiero sacar ccantidad de vecinos etiquetas (de 0 a 9). Sacar el maximo de ahí.
  }

  // El vector y tiene en la i-esima posicion la norma |imagen - y[i]|_2

  std::vector<int> labels(10, 0);

  while (vecinos > 0) {
    int i = indiceMinimo(y); // Me fijo cual es la imagen que minimiza la norma en cada iteracion
    int etiqueta = imagenesTrain.dameEtiqueta(i); // Me fijo cual es la etiqueta de dicho minimo
    labels[etiqueta]++;

    vecinos--;
  }

  return indiceMaximo(labels); // Devuelvo el "más votado"
 
}

int KNN(Matriz& imagenesTrain, Matriz& imagenesTest, int vecinos) {

  int filas = imagenesTest.dimensionFilas();
  int cantidadDeAciertos = 0;
  
  for (int i = 0; i < filas; ++i) {

    int etiqueta = dameEtiqueta(imagenesTrain, imagenesTest[i], vecinos); // Le asigna a qué número pertenece la i-ésima imagen de imagenesTest
    imagenesTest.estimar(i, etiqueta); // En matriz.estimar tengo lo que supongo que es la imagen. En matriz.etiqueta tengo lo que de verdad es

 /*   
    std::cout << "la etiqueta de la imagen " << i << " es " << imagenesTest.dameEtiqueta(i) << std::endl;
    std::cout << "su estimación fue: " << imagenesTest.dameEstimacion(i) << std::endl;

    std::cout << "etiqueta: " << imagenesTest.dameEtiqueta(i) << std::endl;
    std::cout << "estimación: " << imagenesTest.dameEstimacion(i) << std::endl;
  */
    if (etiqueta == imagenesTest.dameEtiqueta(i))
      cantidadDeAciertos++;
  }

  //std::cout << "cantidad de aciertos " << cantidadDeAciertos << std::endl;

  return cantidadDeAciertos; // Esto no sé si es necesario aún

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

//  train = path.append("Test.csv"); 
  train = path.append("train.csv"); 
//  train = path.append("train2.csv"); 
//  train = path.append("train3.csv"); 
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
//    int tamImagen = 1;
//    vector< vector<double> > imagenesTrain2(cantImagenesTrain, vector<double>(tamImagen));
//    vector< vector<double> > imagenesTest2(cantImagenesTest, vector<double>(tamImagen));
//    vector< int > agenesTest2(cantImagenesTrain);
//    vector< int > agenesTest3(cantImagenesTrain);
//    vector< int > agenesTest4(cantImagenesTest);
//    vector< int > agenesTest5(cantImagenesTest);
    Matriz imagenesTrain(cantImagenesTrain, tamImagen);
    Matriz imagenesTest(cantImagenesTest, tamImagen);

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
    // isTrain tiene para cada imagen (de 0 a cantImagenesTotales) si es parte del train o no
    
    switch(method) {
      case 0: { // Método KNN

                
//          imagenesTest.mostrar();
          std::cout << "--------------------------" << std::endl;
//          std::cout << imagenesTrain.dimensionFilas() << std::endl;
//          std::cout << imagenesTrain.dimensionColumnas() << std::endl;
//          std::cout << imagenesTest.dimensionFilas() << std::endl;
//          std::cout << imagenesTest.dimensionColumnas() << std::endl;


//          imagenesTrain.mostrar();
          KNN(imagenesTrain, imagenesTest, vecinos);

          break;
              
        }

      case 1: { // Método KNN+PCA

              
          break;
        }

      case 2: { // Método KNN+PLS-DA

        Matriz& imagenesTrainReducida = PLSDAMethod(imagenesTrain, dimensiones); //Por ahora le hardcodeo el segundo parametro TODO: ver como cambiarlo

        // Hay que preguntar si se hace exactamente lo mismo con los dos o no.
        Matriz& imagenesTestReducida = PLSDAMethod(imagenesTest, dimensiones); //Por ahora le hardcodeo el segundo parametro TODO: ver como cambiarlo
        
        // Para cada imagen en imagenesTests aplicarle la transformación característica (para reducir su dimensión)

        KNN(imagenesTrainReducida, imagenesTestReducida, vecinos);
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
