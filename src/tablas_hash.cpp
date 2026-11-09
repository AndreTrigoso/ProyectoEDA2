#pragma once
#include <iostream>
#include <string>
#include "parser.cpp"

using namespace std;

// ============================================================================
// FUNCIONES MATEMATICAS AUXILIARES (Primos m y Potencias 4^n)
// ============================================================================

bool esPrimo(int n) {
    if (n <= 1) return false;
    if (n == 2) return true;
    if (n % 2 == 0) return false;
    for (int i = 3; i * i <= n; i += 2) {
        if (n % i == 0) return false;
    }
    return true;
}

// Busca el siguiente numero primo m mayor o igual a 'objetivo'
int siguientePrimo(int objetivo) {
    if (objetivo <= 2) return 2;
    int candidato = (objetivo % 2 == 0) ? (objetivo + 1) : objetivo;
    while (!esPrimo(candidato)) {
        candidato += 2;
    }
    return candidato;
}

// Calcula 4^n para el tamano de la Tabla Hash 2 (Double Hashing)
int potenciaDeCuatro(int n) {
    int resultado = 1;
    for (int i = 0; i < n; i++) {
        resultado *= 4;
    }
    return resultado;
}

// ============================================================================
// b.1) TABLA HASH 1: LINEAR PROBING (Almacena ID -> Secuencia de ADN)
// Empieza con table size = 10. Si el load factor excede 1, incrementa al
// siguiente numero primo m que reduzca el load factor a ~0.5.
// ============================================================================

struct CasillaSeq {
    DNASequence dato;
    bool ocupado;
    int sondeos;

    CasillaSeq() {
        ocupado = false;
        sondeos = 0;
    }
};

struct InfoRehash {
    int capAnterior;
    int capNueva;
    int elementos;
    double nuevoLoadFactor;
};

struct TablaHashLineal {
    CasillaSeq* tabla;
    int capacidad;
    int elementos;
    long long colisionesTotales;

    InfoRehash historialRehash[32];
    int cantidadRehashes;

    TablaHashLineal(int tamInicial = 10) {
        capacidad = tamInicial;
        elementos = 0;
        colisionesTotales = 0;
        cantidadRehashes = 0;
        tabla = new CasillaSeq[capacidad];
    }

    ~TablaHashLineal() {
        delete[] tabla;
    }

    // Funcion hash polinomial simple para convertir el string ID en un entero
    unsigned long funcionHash(const string& clave) const {
        unsigned long hash = 0;
        for (int i = 0; i < (int)clave.length(); i++) {
            hash = (hash * 31) + (unsigned char)clave[i];
        }
        return hash % capacidad;
    }

    double obtenerLoadFactor() const {
        return (double)elementos / (double)capacidad;
    }

    // Redimensiona la tabla al primo m cuando el load factor supera 1
    void redimensionar(int nuevaCapacidadPrimo) {
        CasillaSeq* tablaAnterior = tabla;
        int capAnterior = capacidad;

        tabla = new CasillaSeq[nuevaCapacidadPrimo];
        capacidad = nuevaCapacidadPrimo;
        elementos = 0;

        for (int i = 0; i < capAnterior; i++) {
            if (tablaAnterior[i].ocupado) {
                unsigned long posInicial = funcionHash(tablaAnterior[i].dato.id);
                for (int paso = 0; paso < capacidad; paso++) {
                    int idx = (posInicial + paso) % capacidad;
                    if (!tabla[idx].ocupado) {
                        tabla[idx].dato = tablaAnterior[i].dato;
                        tabla[idx].ocupado = true;
                        tabla[idx].sondeos = paso + 1;
                        elementos++;
                        break;
                    }
                    colisionesTotales++;
                }
            }
        }
        delete[] tablaAnterior;
    }

    void insertar(const DNASequence& seq) {
        // Si el nuevo elemento hace que el load factor exceda 1.0:
        // buscamos el siguiente primo m >= 2 * (elementos + 1) para que LF quede ~0.5
        int nuevosElementos = elementos + 1;
        double lfProyectado = (double)nuevosElementos / (double)capacidad;

        if (lfProyectado > 1.0) {
            int capVieja = capacidad;
            int primoM = siguientePrimo(nuevosElementos * 2);
            redimensionar(primoM);

            if (cantidadRehashes < 32) {
                historialRehash[cantidadRehashes].capAnterior = capVieja;
                historialRehash[cantidadRehashes].capNueva = primoM;
                historialRehash[cantidadRehashes].elementos = nuevosElementos;
                historialRehash[cantidadRehashes].nuevoLoadFactor = (double)nuevosElementos / (double)primoM;
                cantidadRehashes++;
            }
        }

        // Linear Probing: h(k, i) = (h(k) + i) % capacidad
        unsigned long posInicial = funcionHash(seq.id);
        for (int i = 0; i < capacidad; i++) {
            int idx = (posInicial + i) % capacidad;
            if (!tabla[idx].ocupado) {
                tabla[idx].dato = seq;
                tabla[idx].ocupado = true;
                tabla[idx].sondeos = i + 1;
                elementos++;
                return;
            }
            if (tabla[idx].dato.id == seq.id) {
                tabla[idx].dato.sequence = seq.sequence;
                return;
            }
            colisionesTotales++;
        }
    }

    const DNASequence* buscar(const string& id, int* sondeosRealizados = nullptr) const {
        unsigned long posInicial = funcionHash(id);
        for (int i = 0; i < capacidad; i++) {
            int idx = (posInicial + i) % capacidad;
            if (!tabla[idx].ocupado) {
                if (sondeosRealizados) *sondeosRealizados = i + 1;
                return nullptr;
            }
            if (tabla[idx].dato.id == id) {
                if (sondeosRealizados) *sondeosRealizados = i + 1;
                return &(tabla[idx].dato);
            }
        }
        if (sondeosRealizados) *sondeosRealizados = capacidad;
        return nullptr;
    }

    void reiniciar() {
        delete[] tabla;
        capacidad = 10;
        elementos = 0;
        colisionesTotales = 0;
        cantidadRehashes = 0;
        tabla = new CasillaSeq[capacidad];
    }
};

// ============================================================================
// b.2) TABLA HASH 2: DOUBLE HASHING (Almacena k-mer -> Frecuencia)
// Tamano redimensionado segun la formula 4^n (ej: 3-mers = 64, 6-mers = 4096).
// ============================================================================

struct KmerConteo {
    string kmer;
    int conteo;

    KmerConteo() {
        kmer = "";
        conteo = 0;
    }
    KmerConteo(const string& k, int c) {
        kmer = k;
        conteo = c;
    }
};

struct CasillaKmer {
    string kmer;
    int conteo;
    bool ocupado;
    int sondeos;

    CasillaKmer() {
        kmer = "";
        conteo = 0;
        ocupado = false;
        sondeos = 0;
    }
};

struct TablaHashDoble {
    CasillaKmer* tabla;
    int longitudK;
    int capacidad; // 4^n bins disponibles
    int kmersUnicos;
    int totalVentanasLeidas;
    int ventanasIgnoradas;
    long long colisionesTotales;

    TablaHashDoble(int k = 3) {
        tabla = nullptr;
        configurarTamano(k);
    }

    ~TablaHashDoble() {
        delete[] tabla;
    }

    // Redimensiona la tabla aplicando la formula 4^n
    void configurarTamano(int nuevoK) {
        if (nuevoK < 1) nuevoK = 1;
        if (nuevoK > 10) nuevoK = 10;
        longitudK = nuevoK;
        capacidad = potenciaDeCuatro(longitudK);
        kmersUnicos = 0;
        totalVentanasLeidas = 0;
        ventanasIgnoradas = 0;
        colisionesTotales = 0;

        delete[] tabla;
        tabla = new CasillaKmer[capacidad];
    }

    // Primera funcion hash: h1(k)
    unsigned long hash1(const string& kmer) const {
        unsigned long h = 0;
        for (int i = 0; i < (int)kmer.length(); i++) {
            h = (h * 37) + (unsigned char)kmer[i];
        }
        return h % capacidad;
    }

    // Segunda funcion hash: h2(k)
    // Retorna siempre un numero impar para ser coprimo con capacidad = 4^n
    unsigned long hash2(const string& kmer) const {
        if (capacidad <= 1) return 1;
        unsigned long h = 0;
        for (int i = 0; i < (int)kmer.length(); i++) {
            h = (h * 53) + (unsigned char)kmer[i] * (i + 1);
        }
        unsigned long pasoImpar = ((h % (capacidad / 2)) * 2) + 1;
        return pasoImpar;
    }

    double obtenerLoadFactor() const {
        return (double)kmersUnicos / (double)capacidad;
    }

    // Inserta un k-mer con Double Hashing: h(k, i) = (h1(k) + i * h2(k)) % 4^n
    void insertarKmer(const string& kmer) {
        if ((int)kmer.length() != longitudK || !esKmerValido(kmer)) {
            ventanasIgnoradas++;
            return;
        }

        unsigned long h1 = hash1(kmer);
        unsigned long h2 = hash2(kmer);

        for (int i = 0; i < capacidad; i++) {
            int idx = (h1 + i * h2) % capacidad;

            if (!tabla[idx].ocupado) {
                tabla[idx].kmer = kmer;
                tabla[idx].conteo = 1;
                tabla[idx].ocupado = true;
                tabla[idx].sondeos = i + 1;
                kmersUnicos++;
                totalVentanasLeidas++;
                return;
            }

            if (tabla[idx].kmer == kmer) {
                tabla[idx].conteo++;
                totalVentanasLeidas++;
                return;
            }

            colisionesTotales++;
        }
    }

    // Lee la secuencia con un desplazamiento de 1 caracter y guarda cada k-mer
    void procesarSecuencia(const string& secuencia) {
        configurarTamano(longitudK);
        int largo = (int)secuencia.length();
        if (largo < longitudK) return;

        for (int pos = 0; pos <= largo - longitudK; pos++) {
            string subcadena = secuencia.substr(pos, longitudK);
            insertarKmer(subcadena);
        }
    }

    // Busca un k-mer usando Double Hashing y devuelve su cantidad de apariciones
    int buscar(const string& kmer, int* sondeosRealizados = nullptr) const {
        if ((int)kmer.length() != longitudK || !esKmerValido(kmer)) {
            if (sondeosRealizados) *sondeosRealizados = 0;
            return 0;
        }

        unsigned long h1 = hash1(kmer);
        unsigned long h2 = hash2(kmer);

        for (int i = 0; i < capacidad; i++) {
            int idx = (h1 + i * h2) % capacidad;
            if (!tabla[idx].ocupado) {
                if (sondeosRealizados) *sondeosRealizados = i + 1;
                return 0;
            }
            if (tabla[idx].kmer == kmer) {
                if (sondeosRealizados) *sondeosRealizados = i + 1;
                return tabla[idx].conteo;
            }
        }
        if (sondeosRealizados) *sondeosRealizados = capacidad;
        return 0;
    }
};
