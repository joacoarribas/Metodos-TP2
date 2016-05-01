#ifndef PLSDA_H_
#define PLDSA_H_

using namespace std;

class PLSDA {

  public:
    PLSDA();
    Matriz& transformacionCaracteristica(Matriz m, vector< vector<double> > w, int n, int dimensiones);
    Matriz& PLSDAMethod(Matriz valores, std::vector<unsigned char> &etiquetas, int dimensiones);
    double multiplicarVectores(vector<double> a, vector<double> b);
    void cargarVector(vector<double> &x);
    void normalizar(vector<double> &x);
};

#endif