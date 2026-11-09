#pragma once
#include <iostream>
#include <string>
#include "parser.cpp"

using namespace std;

// ============================================================================
// FUNCIONES MATEMATICAS BASICAS
// ============================================================================

bool esPrimo(int n) {
    if (n <= 1) return false;
    for (int i = 2; i * i <= n; i++) {
        if (n % i == 0) return false;
    }
    return true;
}

int siguientePrimo(int numero) {
    while (!esPrimo(numero)) {
        numero++;
    }
    return numero;
}

int potenciaDeCuatro(int n) {
    int resultado = 1;
    for (int i = 0; i < n; i++) {
        resultado *= 4;
    }
    return resultado;
}

// ============================================================================
// b.1) TABLA HASH 1: LINEAR PROBING (ID -> Secuencia)
// ============================================================================

struct CasillaSeq {
    DNASequence dato;
    bool ocupado;

    CasillaSeq() {
        ocupado = false;
    }
};

struct TablaHashLineal {
    CasillaSeq* tabla;
    int capacidad;
    int elementos;

    TablaHashLineal(int tamInicial = 10) {
        capacidad = tamInicial;
        elementos = 0;
        tabla = new CasillaSeq[capacidad];
    }

    ~TablaHashLineal() {
        delete[] tabla;
    }

    unsigned long funcionHash(const string& id) const {
        unsigned long hash = 0;
        for (int i = 0; i < (int)id.length(); i++) {
            hash = (hash * 31) + id[i];
        }
        return hash % capacidad;
    }

    double obtenerLoadFactor() const {
        return (double)elementos / (double)capacidad;
    }

    // Crea una tabla del nuevo tamano primo m y reinserta las secuencias
    void redimensionar(int nuevoTamanoPrimo) {
        CasillaSeq* tablaVieja = tabla;
        int capVieja = capacidad;

        capacidad = nuevoTamanoPrimo;
        elementos = 0;
        tabla = new CasillaSeq[capacidad];

        for (int i = 0; i < capVieja; i++) {
            if (tablaVieja[i].ocupado) {
                insertar(tablaVieja[i].dato);
            }
        }
        delete[] tablaVieja;
    }

    void insertar(const DNASequence& seq) {
        // Si al insertar el load factor supera 1, redimensionar al siguiente
        // primo m que baje el load factor a ~0.5 (m = siguientePrimo(2 * nuevoTotal))
        if ((double)(elementos + 1) / capacidad > 1.0) {
            int primoM = siguientePrimo((elementos + 1) * 2);
            redimensionar(primoM);
        }

        unsigned long posInicial = funcionHash(seq.id);

        // Linear Probing: avanza de 1 en 1 si la casilla esta ocupada
        for (int i = 0; i < capacidad; i++) {
            int idx = (posInicial + i) % capacidad;

            if (!tabla[idx].ocupado) {
                tabla[idx].dato = seq;
                tabla[idx].ocupado = true;
                elementos++;
                return;
            }
            if (tabla[idx].dato.id == seq.id) {
                tabla[idx].dato.sequence = seq.sequence;
                return;
            }
        }
    }

    const DNASequence* buscar(const string& id) const {
        unsigned long posInicial = funcionHash(id);

        for (int i = 0; i < capacidad; i++) {
            int idx = (posInicial + i) % capacidad;

            if (!tabla[idx].ocupado) {
                return nullptr; // Casilla vacia: no existe
            }
            if (tabla[idx].dato.id == id) {
                return &(tabla[idx].dato);
            }
        }
        return nullptr;
    }
};

// ============================================================================
// b.2) TABLA HASH 2: DOUBLE HASHING (k-mer -> Frecuencia, tamano 4^n)
// ============================================================================

struct CasillaKmer {
    string kmer;
    int conteo;
    bool ocupado;

    CasillaKmer() {
        kmer = "";
        conteo = 0;
        ocupado = false;
    }
};

struct TablaHashDoble {
    CasillaKmer* tabla;
    int longitudK;
    int capacidad;
    int kmersUnicos;

    TablaHashDoble(int k = 3) {
        tabla = nullptr;
        configurarTamano(k);
    }

    ~TablaHashDoble() {
        delete[] tabla;
    }

    // Redimensiona la tabla segun la formula 4^n
    void configurarTamano(int nuevoK) {
        longitudK = nuevoK;
        capacidad = potenciaDeCuatro(longitudK);
        kmersUnicos = 0;

        delete[] tabla;
        tabla = new CasillaKmer[capacidad];
    }

    // Primera funcion hash: posicion inicial
    unsigned long hash1(const string& kmer) const {
        unsigned long h = 0;
        for (int i = 0; i < (int)kmer.length(); i++) {
            h = (h * 37) + kmer[i];
        }
        return h % capacidad;
    }

    // Segunda funcion hash: tamano del salto (siempre impar para recorrer 4^n)
    unsigned long hash2(const string& kmer) const {
        unsigned long h = 0;
        for (int i = 0; i < (int)kmer.length(); i++) {
            h = (h * 53) + kmer[i];
        }
        return ((h % (capacidad / 2)) * 2) + 1;
    }

    void insertarKmer(const string& kmer) {
        if ((int)kmer.length() != longitudK || !esKmerValido(kmer)) return;

        unsigned long h1 = hash1(kmer);
        unsigned long h2 = hash2(kmer);

        // Double Hashing: (h1 + i * h2) % capacidad
        for (int i = 0; i < capacidad; i++) {
            int idx = (h1 + i * h2) % capacidad;

            if (!tabla[idx].ocupado) {
                tabla[idx].kmer = kmer;
                tabla[idx].conteo = 1;
                tabla[idx].ocupado = true;
                kmersUnicos++;
                return;
            }
            if (tabla[idx].kmer == kmer) {
                tabla[idx].conteo++; // Si ya existe, suma 1 a su frecuencia
                return;
            }
        }
    }

    // Ventana deslizante de 1 caracter sobre la secuencia
    void procesarSecuencia(const string& secuencia) {
        configurarTamano(longitudK);
        int largo = (int)secuencia.length();

        for (int pos = 0; pos <= largo - longitudK; pos++) {
            string subcadena = secuencia.substr(pos, longitudK);
            insertarKmer(subcadena);
        }
    }

    int buscar(const string& kmer) const {
        if ((int)kmer.length() != longitudK || !esKmerValido(kmer)) return 0;

        unsigned long h1 = hash1(kmer);
        unsigned long h2 = hash2(kmer);

        for (int i = 0; i < capacidad; i++) {
            int idx = (h1 + i * h2) % capacidad;

            if (!tabla[idx].ocupado) return 0;
            if (tabla[idx].kmer == kmer) return tabla[idx].conteo;
        }
        return 0;
    }
};
