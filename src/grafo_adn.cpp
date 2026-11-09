#pragma once
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

// ============================================================================
// d) GRAFO DE SECUENCIAS DE ADN (Khan et al., 2024)
// Divide la secuencia en pares de nucleotidos (2-mers) y calcula la matriz de
// adyacencia ponderada con la formula: d = 1 / (j - i)^beta
// ============================================================================

struct ParPosicion {
    int posicion; // Posicion 1, 2, 3, ...
    string par;   // 2-mer: ej. "TC", "AT", "TA", ...
    int idCanonico;
};

struct GrafoADN {
    ParPosicion* listaPares;
    int cantidadPares;
    double beta;

    // Vertices que aparecen en la secuencia (en orden de aparicion)
    string verticesActivos[16];
    int canonicoAActivo[16];
    int cantidadVerticesActivos;

    // Matriz de adyacencia de los vertices presentes en la secuencia
    double matrizActiva[16][16];

    // Matriz de adyacencia completa de los 16 posibles 2-mers (AA, AC, ..., TT)
    double matrizCanonica[16][16];

    GrafoADN() {
        listaPares = nullptr;
        cantidadPares = 0;
        beta = 1.0;
        cantidadVerticesActivos = 0;
        limpiar();
    }

    ~GrafoADN() {
        delete[] listaPares;
    }

    void limpiar() {
        delete[] listaPares;
        listaPares = nullptr;
        cantidadPares = 0;
        cantidadVerticesActivos = 0;

        for (int i = 0; i < 16; i++) {
            verticesActivos[i] = "";
            canonicoAActivo[i] = -1;
            for (int j = 0; j < 16; j++) {
                matrizActiva[i][j] = 0.0;
                matrizCanonica[i][j] = 0.0;
            }
        }
    }

    // Convierte A->0, C->1, G->2, T->3
    static int indiceBase(char c) {
        if (c == 'A') return 0;
        if (c == 'C') return 1;
        if (c == 'G') return 2;
        if (c == 'T') return 3;
        return -1;
    }

    // Convierte un 2-mer ("AA".."TT") en su indice de 0 a 15
    static int parAIndiceCanonico(const string& par) {
        if (par.length() != 2) return -1;
        int b1 = indiceBase(par[0]);
        int b2 = indiceBase(par[1]);
        if (b1 < 0 || b2 < 0) return -1;
        return b1 * 4 + b2;
    }

    // Convierte un indice de 0 a 15 en su string ("AA".."TT")
    static string indiceCanonicoAPar(int idx) {
        const char BASES[4] = {'A', 'C', 'G', 'T'};
        if (idx < 0 || idx >= 16) return "??";
        string s = "  ";
        s[0] = BASES[idx / 4];
        s[1] = BASES[idx % 4];
        return s;
    }

    void construirGrafo(const string& secuencia, double valorBeta = 1.0) {
        limpiar();
        beta = valorBeta;

        int maxPares = (int)secuencia.length() / 2;
        if (maxPares <= 0) return;

        listaPares = new ParPosicion[maxPares];
        cantidadPares = 0;

        // Paso 2 del algoritmo: Dividir la secuencia S en pares de nucleotidos (2-mers)
        for (int p = 0; p < maxPares; p++) {
            string par = secuencia.substr(p * 2, 2);
            if (par[0] >= 'a' && par[0] <= 'z') par[0] = par[0] - 'a' + 'A';
            if (par[1] >= 'a' && par[1] <= 'z') par[1] = par[1] - 'a' + 'A';

            int idCan = parAIndiceCanonico(par);
            if (idCan < 0) continue; // Ignorar par si tiene caracter ambiguo

            listaPares[cantidadPares].posicion = p + 1; // Posicion 1, 2, 3...
            listaPares[cantidadPares].par = par;
            listaPares[cantidadPares].idCanonico = idCan;
            cantidadPares++;

            // Paso 3: Si el vertice no existe aun en G, agregarlo
            if (canonicoAActivo[idCan] == -1) {
                int nuevoIdx = cantidadVerticesActivos++;
                verticesActivos[nuevoIdx] = par;
                canonicoAActivo[idCan] = nuevoIdx;
            }
        }

        // Pasos 4 y 5: Determinar peso d = 1 / (j - i)^beta para cada par (u, v)
        // y construir la matriz de adyacencia ponderada A
        for (int i = 0; i < cantidadPares; i++) {
            int posI = listaPares[i].posicion;
            int uCan = listaPares[i].idCanonico;
            int uAct = canonicoAActivo[uCan];

            for (int j = i + 1; j < cantidadPares; j++) {
                int posJ = listaPares[j].posicion;
                int vCan = listaPares[j].idCanonico;
                int vAct = canonicoAActivo[vCan];

                double distancia = (double)(posJ - posI);
                double peso = 1.0 / pow(distancia, beta);

                matrizActiva[uAct][vAct] += peso;
                matrizCanonica[uCan][vCan] += peso;
            }
        }
    }

    double obtenerPeso(const string& origen, const string& destino) const {
        int u = parAIndiceCanonico(origen);
        int v = parAIndiceCanonico(destino);
        if (u < 0 || v < 0) return 0.0;
        return matrizCanonica[u][v];
    }
};
