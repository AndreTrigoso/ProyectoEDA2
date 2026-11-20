#include "../include/grafo_adn.hpp"
#include <cmath>
#include <cctype>


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
    if (numVertices >= 16) return -1;
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
    int* vertIdx = new int[numPares];
    for (int i = 0; i < numPares; i++) {
        string par = secuencia.substr(i * 2, 2);
        for (int c = 0; c < (int)par.length(); c++) {
            par[c] = static_cast<char>(toupper(static_cast<unsigned char>(par[c])));
        }
        if (!esKmerValido(par)) continue;

        int v = obtenerOAgregarVertice(par);
        if (v < 0) continue;

        listaPares[cantidadPares].posicion = i + 1;
        listaPares[cantidadPares].par = par;
        vertIdx[cantidadPares] = v;
        cantidadPares++;
    }

    // Calcular pesos d = 1 / (j - i)^beta para cada pareja i < j
    for (int i = 0; i < cantidadPares; i++) {
        int u = vertIdx[i];
        int posI = listaPares[i].posicion;

        for (int j = i + 1; j < cantidadPares; j++) {
            int v = vertIdx[j];
            int posJ = listaPares[j].posicion;

            double dist = (double)(posJ - posI);
            double w = 1.0 / pow(dist, beta);
            matriz[u][v] += w;
        }
    }
    delete[] vertIdx;
}

double GrafoADN::obtenerPeso(const string& origen, const string& destino) const {
    string origMayus = origen;
    for (int c = 0; c < (int)origMayus.length(); c++) {
        origMayus[c] = static_cast<char>(toupper(static_cast<unsigned char>(origMayus[c])));
    }
    string destMayus = destino;
    for (int c = 0; c < (int)destMayus.length(); c++) {
        destMayus[c] = static_cast<char>(toupper(static_cast<unsigned char>(destMayus[c])));
    }

    int u = -1;
    int v = -1;

    for (int k = 0; k < numVertices; k++) {
        if (vertices[k] == origMayus) u = k;
        if (vertices[k] == destMayus) v = k;
    }

    if (u == -1 || v == -1) return 0.0;
    return matriz[u][v];
}
