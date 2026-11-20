#ifndef MAX_HEAP_HPP
#define MAX_HEAP_HPP

#include <string>
#include "tablas_hash.hpp"

using namespace std;

// ============================================================================
// c) MAX-HEAP PARA RETORNAR LOS N PATRONES MAS FRECUENTES
// ============================================================================

struct NodoHeap {
    string kmer;
    int conteo;
};

struct MaxHeapKmers {
    NodoHeap* heap;
    int tamano;
    int capacidad;

    // ------------------------------------------------------------------------
    // [GENERADO POR IA PARA CONTRASTACION EN LA DISCUSION]
    // Prompt utilizado: "Agrega un contador comparacionesHeap en MaxHeapKmers
    // para medir las comparaciones realizadas al construir y extraer el Top-N
    // y contrastarlo frente a un ordenamiento total en la discusion."
    // ------------------------------------------------------------------------
    int comparacionesHeap;

    MaxHeapKmers();
    ~MaxHeapKmers();
    MaxHeapKmers(const MaxHeapKmers&) = delete;
    MaxHeapKmers& operator=(const MaxHeapKmers&) = delete;

    void intercambiar(int i, int j);
    void subir(int idx);
    void bajar(int idx);
    void insertar(const string& kmer, int conteo);
    void construirDesdeTablaHash(const TablaHashDoble& tablaHash2);
    NodoHeap extraerMaximo();
};

#endif
