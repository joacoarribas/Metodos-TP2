#include "clases/Matriz.h"
//#include "metodos/PLSDA.cpp"
//#include "metodos/PLSDATest.cpp"
//#include "metodos/pca.cpp"
#include "metodos/metodoPotencia.cpp"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream> 
#include <algorithm>
#include <map>
#include <queue>
#include <sys/time.h>
#include <limits>

/* Empieza funciones Toma de tiempos */

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

/* Termina funciones Toma de tiempos */

/* Empieza Auxiliares estadística */

void calcularEstadisticas(Matriz& imagenes, std::ofstream& fEstadisticas, std::ofstream& fEstPrecision, std::ofstream& fEstRecall, std::ofstream& fEstF1) {
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

  for (int i = 0; i < 10; ++i) {
    double precision = (double)tp[i] / (double)(tp[i] + fp[i]);
    precisiones[i] = precision;

    double recall = (double)tp[i] / (double)(tp[i] + fn[i]);
    recalls[i] = recall;
  }

  fEstadisticas << "-----------------------------------------" << std::endl;
  fEstadisticas << "Cálculo precision: " << std::endl;

  double precision = 0;
  for (int i = 0; i < 10; ++i) {
    fEstadisticas << "Precisión clase " << i << ": " << std::fixed << precisiones[i] << std::endl; 
    fEstPrecision << std::fixed << precisiones[i] << std::endl;
    precision += precisiones[i];
  }

  precision /= 10.0;

  fEstadisticas << "Precisión categorizador: " << std::fixed << precision << std::endl; 
  fEstPrecision << std::fixed << precision << std::endl;

  fEstadisticas << "-----------------------------------------" << std::endl;
  fEstadisticas << "Cálculo recall: " << std::endl;

  double recall = 0;
  for (int i = 0; i < 10; ++i) {
    fEstadisticas << "Recall clase " << i << ": " << std::fixed << recalls[i] << std::endl; 
    fEstRecall << std::fixed << recalls[i] << std::endl;
    recall += recalls[i];
  }

  recall /= 10.0;

  fEstadisticas << "Recall categorizador: " << std::fixed << recall << std::endl; 
  fEstRecall << std::fixed << recall << std::endl;

  double F1 = (2 * precision * recall) / (precision + recall);

  fEstadisticas << "F1: " << std::fixed << F1 << std::endl; 
  fEstF1 << std::fixed << F1 << std::endl;
  fEstadisticas << "-----------------------------------------" << std::endl;
  fEstadisticas << "-----------------------------------------" << std::endl;

}

/* Termina auxiliares Estadistica */

/* Empieza auxiliares metodos */

double dameNorma(std::vector<double>& x) {
  double norm = 0;
  int length = x.size();
  for (int i = 0; i < length; ++i)
    norm += x[i] * x[i];
    
  return sqrt(norm);
}

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

void calcular_media(Matriz& matriz, vector<double>& v) { //ver si se pierde precisión con la división
  int n = matriz.dimensionFilas();
  int m = matriz.dimensionColumnas();

  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i)
      v[j] += matriz[i][j]; //acumulo

    v[j] /= (double)n; //calculo promedio final
  }
}

void calcular_matriz_X(Matriz& matriz, Matriz& res){
  int n = matriz.dimensionFilas();
  int m = matriz.dimensionColumnas();

  vector<double> media(m, 0);
  
  calcular_media(matriz, media);
  
  double divisor = sqrt ((double)n - 1.0);

  for (int j = 0; j < m; ++j)
    for (int i = 0; i < n; ++i)
      res[i][j] = (matriz[i][j] - media[j]) / divisor; //escribo resultado en la matriz

}

/* Termina auxiliares Metodos */

/* Empieza PCA */

void calcular_base_ortonormal(Matriz& matriz, Matriz& matriz_ortonormal, int alfa, std::ofstream& filewrite){ //deja en matriz_ortonormal una matriz de alfa columnas
  int m = matriz_ortonormal.dimensionColumnas();

  for (int i = 0; i < alfa; ++i) { //repito alfa veces (hay que experimentar con dicho valor)

    vector<double>& aux = matriz_ortonormal[i];

    Matriz::cargarVector(aux);

    double autovalor = metodoPotencia(matriz, aux); //calculo el i-ésimo autovalor, en aux queda el autovector
    filewrite << std::scientific << autovalor << std::endl;

    normalizar(aux);

    /* Deflación */
    Matriz auxiliar(m, m);

    vector<double> aux2(m, 0);

    for (int i = 0; i < m; ++i) {
      aux2[i] = aux[i] * autovalor;
    }

    auxiliar.multiplicarVectoresDameMatriz(aux, aux2);
    matriz.menos(auxiliar);

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

  Xt.multiplicarMatrices(X, m_covarianza); //crear matriz covarianza

  /* Reducción de la dimensión */
  calcular_base_ortonormal(m_covarianza, m_ortonormal, alfa, filewrite); //deja en matriz_ortonormal una matriz de alfa filas

  /* Transformación característica */
  Matriz m_ortonormal_traspuesta(m, alfa);
  m_ortonormal.trasponer(m_ortonormal_traspuesta);

  /* Transformación característica */
  matriz.multiplicarMatrices(m_ortonormal_traspuesta, res); // Chequear esto en caso de que no ande
}

/* Termina PCA */


/* Empieza PLSDA */

void PLSDAMethod(Matriz& imagenes, Matriz& imagenesTrainPLSDA, Matriz& autovectores, int dimensiones, std::ostream& filewrite) {

	long n = imagenes.dimensionFilas();
	long m = imagenes.dimensionColumnas();
	double raizN = sqrt((double)n-1.0);

	Matriz X(n, m);
	Matriz Y(n, 10);

	//calcular mu = promedio de la suma de todos los vectores
  std::vector<double> promedios(n, 0);
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i)
			promedios[j] += imagenes[i][j];

		promedios[j] /= (double)n; 
	}

	//X matriz de nXm donde cada fila es traspuesta(x_i − μ)/raiz(n − 1)
  for (int j = 0; j < m; ++j)
    for (int i = 0; i < n; ++i)
			X[i][j] = (imagenes[i][j] - promedios[j]) / raizN;

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
  std::vector<double> promY(10, 0);
  for (int j = 0; j < 10; ++j) {
    for (int i = 0; i < n; ++i)
			promY[j] += preY[i][j];

		promY[j] /= (double)n; //esto puede traer error
	}

	//defino Y como preY_i − promY_i / raiz(n-1)
  for (int j = 0; j < 10; ++j)
    for (int i = 0; i < n; ++i)
			Y[i][j] = (preY[i][j] - promY[j]) / raizN;

  // Tengo dimensiones cantidad de autovectores de tamaño m

	for (int i=0; i<dimensiones; ++i) {

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

		//calcular el autovector asociado al mayor autovalo de Mi
		std::vector<double>& autovector = autovectores[i]; // Si w = nxm entonces wi = m
		Matriz::cargarVector(autovector);

    double autovalor = metodoPotencia(Mi, autovector); //descarto el autovalor que vino, solo necesito el auvector en wi
    filewrite << std::scientific << autovalor << std::endl;

		//normalizar autovector con norma 2
		normalizar(autovector);
		
		//ti = X * wi
		std::vector<double> ti(n, 0);
    X.multiplicarVectorDer(autovector, ti); // x = nxm y wi = mx1 entonces ti = nx1; 
		
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
    std::cout << std::endl;
	}

  // Los autovectores estan como FILAS de la matriz w. Para multiplicarla por imagenes tengo que trasponerla
  
  Matriz wTras(m, dimensiones);
  autovectores.trasponer(wTras);

  imagenes.multiplicarMatrices(wTras, imagenesTrainPLSDA); // Esto es la transformación caracteritisca

}

/* Termina PLSDA */


/* Empieza KNN */
    
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

int KNN(Matriz& imagenesTrain, Matriz& imagenesTest, int vecinos,
    std::ofstream& fEstadisticas, std::ofstream& fEstPresicion,
    std::ofstream& fEstRecall, std::ofstream& fEstF1, std::ofstream& fEstHitRate) {

  std::cout << "entro knn" << std::endl;

  int filas = imagenesTest.dimensionFilas();
  int hitRate = 0;

  for (int i = 0; i < filas; ++i) {

    int etiqueta = dameEtiquetaEstimada(imagenesTrain, imagenesTest[i], vecinos); // Le asigna a qué número pertenece la i-ésima imagen de imagenesTest
    imagenesTest.estimar(i, etiqueta); // En matriz.estimar tengo lo que supongo que es la imagen. En matriz.etiqueta tengo lo que de verdad es

    if (etiqueta == imagenesTest.dameEtiqueta(i))
      hitRate++;
  }

  std::cout << "termino knn" << std::endl;
  fEstadisticas << "cantidad de imagenes Test: " << filas << std::endl;
  fEstadisticas << "cantidad de aciertos (Hit Rate): " << hitRate << std::endl;
  fEstHitRate << hitRate << std::endl;
  acum += get_time();
  fEstadisticas << "tiempo: " << acum << std::endl;
  acum = 0;

  calcularEstadisticas(imagenesTest, fEstadisticas, fEstPresicion, fEstRecall, fEstF1);

  std::cout << "Hit rate: " << hitRate << std::endl;
  return hitRate; // Esto no sé si es necesario aún

}

/* Termina Knn */

void evaluarTests(std::string fileTestData, std::ofstream& fileTestResult, std::ofstream& fileEstadisticas, std::ofstream& fEstPrecision, std::ofstream& fEstRecall, std::ofstream& fEstF1, std::ofstream& fEstHitRate, std::string path, int method, int vecinos, int componentes, int dimensiones, int particiones) {

  std::string lineData;
  std::string lineTest;
  std::string lineTrain;
  std::ifstream fileData (fileTestData.c_str());
  int z = 0;

  getline(fileData, lineData); // Pido primer línea para instanciar todo
  std::istringstream issData(lineData);

  std::string train;
  std::string test;

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
    Matriz autovectoresPLSDA(dimensiones, tamImagen);
    init_time();
 
    if (method == 1) {
	    std::cout << "entre aca | componentes -> " << componentes << std::endl; 
      PCAMethod(imagenesTrain, imagenesTrainPCA, autovectoresPCA, componentes, fileTestResult);
      std::cout << "sali ah reee" << std::endl;
    } else {
      if (method == 2) {
	std::cout << "entre aca | dimensiones -> " << dimensiones << std::endl; 
      PLSDAMethod(imagenesTrain, imagenesTrainPLSDA, autovectoresPLSDA, dimensiones, fileTestResult); // Por ahora le hardcodeo el segundo parametro TODO: ver como cambiarloi
	std::cout << "sali ah reee" << std::endl;
      }
    }
    
    fileEstadisticas << "Iteración: " << z << std::endl;
    fileEstadisticas << "Vecinos: " << vecinos << std::endl;
    fileEstadisticas << "dimensiones: " << dimensiones << std::endl;
    fileEstadisticas << "componentes: " << componentes << std::endl;
    fileEstadisticas << "particiones: " << particiones << std::endl;

    switch(method) {
      case 0: { // Método KNN

          KNN(imagenesTrain, imagenesTest, vecinos, fileEstadisticas, fEstPrecision, fEstRecall, fEstF1, fEstHitRate);

          break;
              
        }

      case 1: { // Método KNN+PCA

          Matriz wTras(tamImagen, componentes);
          autovectoresPCA.trasponer(wTras);

          imagenesTest.multiplicarMatrices(wTras, imagenesTestPCA);

          KNN(imagenesTrainPCA, imagenesTestPCA, vecinos, fileEstadisticas, fEstPrecision, fEstRecall, fEstF1, fEstHitRate);
                
          break;
        }

      case 2: { // Método KNN+PLS-DA

        Matriz wTras(tamImagen, dimensiones);
        autovectoresPLSDA.trasponer(wTras);

        imagenesTest.multiplicarMatrices(wTras, imagenesTestPLSDA);

        KNN(imagenesTrainPLSDA, imagenesTestPLSDA, vecinos, fileEstadisticas, fEstPrecision, fEstRecall, fEstF1, fEstHitRate);
        break;

        }

      case 3: { // Método KNN+PLS-DA

        Matriz wTras(tamImagen, dimensiones);
        autovectoresPLSDA.trasponer(wTras);

        imagenesTest.multiplicarMatrices(wTras, imagenesTestPLSDA);

        KNN(imagenesTrainPLSDA, imagenesTestPLSDA, vecinos, fileEstadisticas, fEstPrecision, fEstRecall, fEstF1, fEstHitRate);
        break;

        }

    }

    ++z;
    std::cout << z << std::endl; // SOlo uso esto para ver cuantas iteraciones de lineas de archivo hizo

    fileTestResult << "ImageId,Label" << std::endl;
    for (int k = 1; k <= cantImagenesTest; ++k) {
      fileTestResult << k << "," << imagenesTestPLSDA.dameEstimacion(k-1) << std::endl;
    }

  }

  acum += get_time();
  std::cout << std::fixed << acum << std::endl;
}

void evaluarTests2(std::string fileTestData, std::ofstream& fileTestResult, std::ofstream& fileEstadisticas, std::ofstream& fEstPrecision, std::ofstream& fEstRecall, std::ofstream& fEstF1, std::ofstream& fEstHitRate, std::string path, int method, int vecinos, int componentes, int dimensiones, int particiones) {

  std::string lineTest;
  std::string lineTrain;
  int z = 0;

  std::string train;
  std::string test;

  train = path.append("train.csv"); 
  test = path.append("test.csv"); 

  std::ifstream fileTrain (train.c_str());
  std::ifstream fileTest (test.c_str()); // Hasta aca sólo instancie variables

  int tamImagen = 784;
  int cantImagenesTest = 28000;
  int cantImagenesTrain = 42000;

  Matriz imagenesTrain(cantImagenesTrain, tamImagen);
  Matriz imagenesTrainPCA(cantImagenesTrain, componentes);
  Matriz imagenesTrainPLSDA(cantImagenesTrain, dimensiones);

  Matriz imagenesTest(cantImagenesTest, tamImagen);
  Matriz imagenesTestPCA(cantImagenesTest, componentes);
  Matriz imagenesTestPLSDA(cantImagenesTest, dimensiones);

  getline (fileTrain, lineTrain); // Descarto la primer línea que sólo tiene strings "label, pixel0, ..."

  for (int k = 0; k < cantImagenesTrain; ++k) {
  
    getline (fileTrain, lineTrain);

    std::istringstream issTrain(lineTrain);

    int j = 0;
    int a = 0;
    std::string s;

    while (getline(issTrain, s, ',')) {
      if (j == 0) {
        imagenesTrain.etiquetar(k, atoi(s.c_str()));
        imagenesTrainPLSDA.etiquetar(k, atoi(s.c_str()));
        imagenesTrainPCA.etiquetar(k, atoi(s.c_str()));
        ++j;
      } else {
        imagenesTrain[k][a] = std::stod(s);
        ++a;
      }
    }

  }

  getline (fileTest, lineTest); // Descarto la primer línea que sólo tiene strings "label, pixel0, ..."

  for (int k = 0; k < cantImagenesTest; ++k) {
  
    getline (fileTest, lineTest);

    std::istringstream issTest(lineTest);

    int j = 0;
    int a = 0;
    std::string s;

    while (getline(issTest, s, ',')) {
      imagenesTest[k][a] = std::stod(s);
      ++a;
    }
  }

    // imagenesTrain tiene en cada fila una imagen para la parte de train Y su correspondiente etiqueta (de 0 a 9)
    // imagenesTest tiene en cada fila una imagen para la parte de test Y su correspondiente etiqueta (de 0 a 9)
    // imagenesTrainPSLA tendrá la transformación característica de PSLA
    // imagenesTrainPCA tendrá la transformación característica de PCA
    // isTrain tiene para cada imagen (de 0 a cantImagenesTotales) si es parte del train o no


    Matriz autovectoresPCA(componentes, tamImagen);
    Matriz autovectoresPLSDA(dimensiones, tamImagen);
    init_time();
 
    if (method == 1) {
	    std::cout << "entre aca | componentes -> " << componentes << std::endl; 
      PCAMethod(imagenesTrain, imagenesTrainPCA, autovectoresPCA, componentes, fileTestResult);
      std::cout << "sali ah reee" << std::endl;
    } else {
      if (method == 2) {
	std::cout << "entre aca | dimensiones -> " << dimensiones << std::endl; 
      PLSDAMethod(imagenesTrain, imagenesTrainPLSDA, autovectoresPLSDA, dimensiones, fileTestResult); // Por ahora le hardcodeo el segundo parametro TODO: ver como cambiarloi
	std::cout << "sali ah reee" << std::endl;
      }
    }
    
    fileEstadisticas << "Iteración: " << z << std::endl;
    fileEstadisticas << "Vecinos: " << vecinos << std::endl;
    fileEstadisticas << "dimensiones: " << dimensiones << std::endl;
    fileEstadisticas << "componentes: " << componentes << std::endl;
    fileEstadisticas << "particiones: " << particiones << std::endl;

    switch(method) {
      case 0: { // Método KNN

          KNN(imagenesTrain, imagenesTest, vecinos, fileEstadisticas, fEstPrecision, fEstRecall, fEstF1, fEstHitRate);

          break;
              
        }

      case 1: { // Método KNN+PCA

          Matriz wTras(tamImagen, componentes);
          autovectoresPCA.trasponer(wTras);

          imagenesTest.multiplicarMatrices(wTras, imagenesTestPCA);

          KNN(imagenesTrainPCA, imagenesTestPCA, vecinos, fileEstadisticas, fEstPrecision, fEstRecall, fEstF1, fEstHitRate);
                
          break;
        }

      case 2: { // Método KNN+PLS-DA

        Matriz wTras(tamImagen, dimensiones);
        autovectoresPLSDA.trasponer(wTras);

        imagenesTest.multiplicarMatrices(wTras, imagenesTestPLSDA);

        KNN(imagenesTrainPLSDA, imagenesTestPLSDA, vecinos, fileEstadisticas, fEstPrecision, fEstRecall, fEstF1, fEstHitRate);
        break;

        }

    }



    ++z;
    std::cout << z << std::endl; // SOlo uso esto para ver cuantas iteraciones de lineas de archivo hizo

  acum += get_time();
  std::cout << std::fixed << acum << std::endl;

  fileTestResult << "ImageId,Label" << std::endl;
  for (int k = 1; k <= cantImagenesTest; ++k) {
    fileTestResult << k << "," << imagenesTestPLSDA.dameEstimacion(k-1) << std::endl;
  }

}



int elegirOpciones(std::string fileTestData, std::string fileTestResult, std::string fileEstadisticas, std::string filePrecision, std::string fileRecall, std::string fileF1, std::string fileHitRate, int method) {
  
  std::string lineData;
  std::string lineTest;
  std::string lineTrain;
  std::ifstream fileData (fileTestData.c_str());
  std::ofstream fileWrite (fileTestResult.c_str());
  std::ofstream fileWrite2 (fileEstadisticas.c_str());
  std::ofstream fileWrite3 (filePrecision.c_str());
  std::ofstream fileWrite4 (fileRecall.c_str());
  std::ofstream fileWrite5 (fileF1.c_str());
  std::ofstream fileWrite6 (fileHitRate.c_str());
  int z = 0;

  getline(fileData, lineData); // Pido primer línea para instanciar todo
  std::istringstream issData(lineData);

  std::string path;
  std::string train;
  std::string test;
  
  std::vector<int> vecinos(3, 0);
  std::vector<int> componentes(3, 0);
  std::vector<int> dimensiones(3, 0);
  std::vector<int> particiones(3, 0);

  int varAModificar;
  std::cout << "¿queres modifcar variables o sar por defecto?" << std::endl;
  std::cout << "1: modificar" << std::endl;
  std::cout << "2: por defecto" << std::endl;
  std::cin >> varAModificar;  
   
  issData >> path;
  
  if (varAModificar == 1) {
    std::cout << "ingrese el rango de vecinos (el intervalo es cerrado)" << std::endl;
    std::cout << "Inicio: ";
    std::cin >> vecinos[0];
    std::cout << "Fin: ";
    std::cin >> vecinos[1];
    std::cout << "Step: ";
    std::cin >> vecinos[2];

    std::cout << "ingrese el rango de las componentes (el intervalo es cerrado)" << std::endl;
    std::cout << "Inicio: ";
    std::cin >> componentes[0];
    std::cout << "Fin: ";
    std::cin >> componentes[1];
    std::cout << "Step: ";
    std::cin >> componentes[2];

    std::cout << "ingrese el rango de las dimensiones (el intervalo es cerrado)" << std::endl;
    std::cout << "Inicio: ";
    std::cin >> dimensiones[0];
    std::cout << "Fin: ";
    std::cin >> dimensiones[1];
    std::cout << "Step: ";
    std::cin >> dimensiones[2];

    std::cout << "ingrese el rango de las particiones para cross folding (el intervalo es cerrado)" << std::endl;
    std::cout << "Inicio: ";
    std::cin >> particiones[0];
    std::cout << "Fin: ";
    std::cin >> particiones[1];
    std::cout << "Step: ";
    std::cin >> particiones[2];

    for (int vec = vecinos[0]; vec <= vecinos[1]; vec+=vecinos[2]) {  
      for (int comp = componentes[0]; comp <= componentes[1]; comp+=componentes[2]) {
        for (int dim = dimensiones[0]; dim <= dimensiones[1]; dim+=dimensiones[2]) {
          for (int part = particiones[0]; part <= particiones[1]; part+=particiones[2]) {
            //evaluarTests(fileTestData, fileWrite, fileWrite2, fileWrite3, fileWrite4, fileWrite5, fileWrite6, path, method, vec, comp, dim, part);
            evaluarTests2(fileTestData, fileWrite, fileWrite2, fileWrite3, fileWrite4, fileWrite5, fileWrite6, path, method, vec, comp, dim, part);
          }
        }
      }
    }

  } else {
    issData >> vecinos[0];
    issData >> componentes[0];
    issData >> dimensiones[0];
    issData >> particiones[0];
    evaluarTests(fileTestData, fileWrite, fileWrite2, fileWrite3, fileWrite4, fileWrite5, fileWrite6, path, method, vecinos[0], componentes[0], dimensiones[0], particiones[0]);
  }
  return 0;
}

int main(int argc, char** argv) {
  std::string fileTestData(argv[1]);
  std::string fileTestResult(argv[2]);
  std::string fileTestEstadisticas(argv[3]);
  std::string fileTestPrecision(argv[4]);
  std::string fileTestRecall(argv[5]);
  std::string fileTestF1(argv[6]);
  std::string fileTestHitRate(argv[7]);
  int method(atoi(argv[8]));
  // Recibo por parametro tres archivos
  // El primero del cual leo los datos a evaluar
  // El segundo en el cual evaluo si los resultados fueron correctos
  // El tercero el método a realizar (0 KNN, 1 PCC+KNN, 2 PLS-DA+KNN)

  elegirOpciones(fileTestData, fileTestResult, fileTestEstadisticas, fileTestPrecision, fileTestRecall, fileTestF1, fileTestHitRate, method);

  return 0;
}
