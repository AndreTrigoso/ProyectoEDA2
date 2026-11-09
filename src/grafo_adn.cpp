#pragma once
#include <iostream>
#include <string>
#include <cmath>
#include "parser.cpp"

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

    GrafoADN() {
        listaPares = nullptr;
        cantidadPares = 0;
        numVertices = 0;
        beta = 1.0;
    }

    ~GrafoADN() {
        delete[] listaPares;
    }

    // Si el vertice 2-mer no existe en G lo agrega; devuelve su indice (0..numVertices-1)
    int obtenerOAgregarVertice(const string& par) {
        for (int i = 0; i < numVertices; i++) {
            if (vertices[i] == par) return i;
        }
        vertices[numVertices] = par;
        numVertices++;
        return numVertices - 1;
    }

    void construirGrafo(const string& secuencia, double valorBeta = 1.0) {
        delete[] listaPares;
        cantidadPares = 0;
        numVertices = 0;
        beta = valorBeta;

        // Paso 1: Iniciar con el grafo G vacio (matriz en ceros)
        for (int i = 0; i < 16; i++) {
            for (int j = 0; j < 16; j++) {
                matriz[i][j] = 0.0;
            }
        }

        int totalPares = (int)secuencia.length() / 2;
        if (totalPares <= 0) return;
        listaPares = new ParPosicion[totalPares];

        // Paso 2 y 3: Dividir S en pares de nucleotidos y registrar sus vertices en G
        for (int p = 0; p < totalPares; p++) {
            string par = secuencia.substr(p * 2, 2);
            if (!esKmerValido(par)) continue;

            listaPares[cantidadPares].posicion = p + 1;
            listaPares[cantidadPares].par = par;
            cantidadPares++;

            obtenerOAgregarVertice(par);
        }

        // Paso 4 y 5: Calcular peso d = 1 / (j - i)^beta y llenar la matriz A[u, v]
        for (int i = 0; i < cantidadPares; i++) {
            int u = obtenerOAgregarVertice(listaPares[i].par);
            int posI = listaPares[i].posicion;

            for (int j = i + 1; j < cantidadPares; j++) {
                int v = obtenerOAgregarVertice(listaPares[j].par);
                int posJ = listaPares[j].posicion;

                double peso = 1.0 / pow((double)(posJ - posI), beta);
                matriz[u][v] += peso;
            }
        }
    }

    double obtenerPeso(const string& origen, const string& destino) const {
        int u = -1, v = -1;
        for (int i = 0; i < numVertices; i++) {
            if (vertices[i] == origen) u = i;
            if (vertices[i] == destino) v = i;
        }
        if (u == -1 || v == -1) return 0.0;
        return matriz[u][v];
    }
};
