#include "../include/grafo_adn.hpp"
#include <cmath>


// d) GRAFO DIRIGIDO DE 2-MERS Y MATRIZ DE ADYACENCIA (Khan et al. 2024)
// - Divide S en pares de nucleotidos disjuntos (paso 2)
// - Calcula d = 1 / (j - i)^beta y llena la matriz de adyacencia A[u][v]


GrafoADN::GrafoADN() {
    listaPares = nullptr;
    cantidadPares = 0;
    numVertices = 0;
    beta = 1.0;

    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 16; j++) {
            matriz[i][j] = 0.0;
        }
    }
}

GrafoADN::~GrafoADN() {
    delete[] listaPares;
}

int GrafoADN::obtenerOAgregarVertice(const string& par) {
    for (int v = 0; v < numVertices; v++) {
        if (vertices[v] == par) {
            return v;
        }
    }
    vertices[numVertices] = par;
    numVertices++;
    return numVertices - 1;
}

void GrafoADN::construirGrafo(const string& secuencia, double valorBeta) {
    delete[] listaPares;
    listaPares = nullptr;
    cantidadPares = 0;
    numVertices = 0;
    beta = valorBeta;

    // Inicializar Matriz de Adyacencia en 0.0 (Semana 8)
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 16; j++) {
            matriz[i][j] = 0.0;
        }
    }

    int numPares = (int)secuencia.length() / 2;
    if (numPares <= 0) return;

    // Arreglo dinamico con los 2-mers en sus posiciones 1..numPares
    listaPares = new ParPosicion[numPares];
    for (int i = 0; i < numPares; i++) {
        string par = secuencia.substr(i * 2, 2);
        if (!esKmerValido(par)) continue;

        listaPares[cantidadPares].posicion = i + 1;
        listaPares[cantidadPares].par = par;
        cantidadPares++;

        obtenerOAgregarVertice(par);
    }

    // Calcular pesos d = 1 / (j - i)^beta para cada pareja i < j
    for (int i = 0; i < cantidadPares; i++) {
        int u = obtenerOAgregarVertice(listaPares[i].par);
        int posI = listaPares[i].posicion;

        for (int j = i + 1; j < cantidadPares; j++) {
            int v = obtenerOAgregarVertice(listaPares[j].par);
            int posJ = listaPares[j].posicion;

            double dist = (double)(posJ - posI);
            double w = 1.0 / pow(dist, beta);
            matriz[u][v] += w;
        }
    }
}

double GrafoADN::obtenerPeso(const string& origen, const string& destino) const {
    int u = -1;
    int v = -1;

    for (int k = 0; k < numVertices; k++) {
        if (vertices[k] == origen) u = k;
        if (vertices[k] == destino) v = k;
    }

    if (u == -1 || v == -1) return 0.0;
    return matriz[u][v];
}