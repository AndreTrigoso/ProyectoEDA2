#ifndef GRAFO_ADN_HPP
#define GRAFO_ADN_HPP

#include <string>
#include "parser.hpp"

using namespace std;

// ============================================================================
// d) GRAFO DE SECUENCIAS DE ADN (Algoritmo de Khan et al., 2024)
// ============================================================================

struct ParPosicion {
    int posicion; // 1, 2, 3...
    string par;   // 2-mer: "TC", "AT", etc.
};

struct GrafoADN {
    ParPosicion* listaPares;
    int cantidadPares;

    string vertices[16];     // Vertices unicos de 2-mers en el grafo G
    int numVertices;
    double matriz[16][16];   // Matriz de adyacencia ponderada A[u][v]
    double beta;

    GrafoADN();
    ~GrafoADN();
    GrafoADN(const GrafoADN&) = delete;
    GrafoADN& operator=(const GrafoADN&) = delete;

    int obtenerOAgregarVertice(const string& par);
    void construirGrafo(const string& secuencia, double valorBeta = 1.0);
    double obtenerPeso(const string& origen, const string& destino) const;
};

#endif
