#ifndef TABLAS_HASH_HPP
#define TABLAS_HASH_HPP

#include <string>
#include "parser.hpp"

using namespace std;

// Funciones matematicas basicas
bool esPrimo(int n);
int siguientePrimo(int numero);
int potenciaDeCuatro(int n);

// ============================================================================
// b.1) TABLA HASH 1: LINEAR PROBING (ID -> Secuencia)
// ============================================================================

struct CasillaSeq {
    DNASequence dato;
    bool ocupado;

    CasillaSeq();
};

struct TablaHashLineal {
    CasillaSeq* tabla;
    int capacidad;
    int elementos;

    // ------------------------------------------------------------------------
    // [GENERADO POR IA PARA CONTRASTACION EN LA DISCUSION]
    // Prompt utilizado: "Agrega variables de instrumentacion en TablaHashLineal
    // para medir sondas totales, colisiones, maximo de sondas, rehashes y
    // sondas de busqueda, con el fin de contrastar linear probing vs double
    // hashing en la discusion del informe."
    // ------------------------------------------------------------------------
    int totalSondas;
    int colisiones;
    int maxSondas;
    int numRehashes;
    mutable int ultimasSondasBusqueda;
    bool enRehash;

    TablaHashLineal(int tamInicial = 10);
    ~TablaHashLineal();

    unsigned long funcionHash(const string& id) const;
    double obtenerLoadFactor() const;
    void redimensionar(int nuevoTamanoPrimo);
    void insertar(const DNASequence& seq);
    const DNASequence* buscar(const string& id) const;
    void reiniciar();
};

// ============================================================================
// b.2) TABLA HASH 2: DOUBLE HASHING (k-mer -> Frecuencia, tamano 4^n)
// ============================================================================

struct CasillaKmer {
    string kmer;
    int conteo;
    bool ocupado;
    int sondas; // [GENERADO POR IA PARA CONTRASTACION] Sondas al insertar este k-mer

    CasillaKmer();
};

struct TablaHashDoble {
    CasillaKmer* tabla;
    int longitudK;
    int capacidad;
    int kmersUnicos;

    // ------------------------------------------------------------------------
    // [GENERADO POR IA PARA CONTRASTACION EN LA DISCUSION]
    // Prompt utilizado: "Agrega variables y un metodo de contrastacion en
    // TablaHashDoble para medir sondas totales, colisiones y maximo de sondas
    // de Double Hashing frente a Linear Probing usando la misma hash1 y
    // capacidad 4^k."
    // ------------------------------------------------------------------------
    int totalSondas;
    int colisiones;
    int maxSondas;
    int sondasLinearProbing;
    int colisionesLinearProbing;
    int maxSondasLinearProbing;
    mutable int ultimasSondasBusqueda;

    TablaHashDoble(int k = 3);
    ~TablaHashDoble();

    void configurarTamano(int nuevoK);
    unsigned long hash1(const string& kmer) const;
    unsigned long hash2(const string& kmer) const;
    void insertarKmer(const string& kmer);
    void procesarSecuencia(const string& secuencia);
    void contrastarConLinearProbing(const string& secuencia); // [GENERADO POR IA]
    int buscar(const string& kmer) const;
};

#endif
