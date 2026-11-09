#pragma once
#include <iostream>
#include <string>
#include "tablas_hash.cpp"

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

    MaxHeapKmers() {
        capacidad = 64;
        tamano = 0;
        heap = new NodoHeap[capacidad];
    }

    ~MaxHeapKmers() {
        delete[] heap;
    }

    void intercambiar(int i, int j) {
        NodoHeap temp = heap[i];
        heap[i] = heap[j];
        heap[j] = temp;
    }

    // Sube el nodo mientras su conteo sea mayor que el de su padre
    void subir(int idx) {
        while (idx > 0) {
            int padre = (idx - 1) / 2;
            if (heap[idx].conteo > heap[padre].conteo) {
                intercambiar(idx, padre);
                idx = padre;
            } else {
                break;
            }
        }
    }

    // Baja el nodo intercambiandolo con el mayor de sus dos hijos
    void bajar(int idx) {
        while (true) {
            int izq = 2 * idx + 1;
            int der = 2 * idx + 2;
            int mayor = idx;

            if (izq < tamano && heap[izq].conteo > heap[mayor].conteo) {
                mayor = izq;
            }
            if (der < tamano && heap[der].conteo > heap[mayor].conteo) {
                mayor = der;
            }

            if (mayor != idx) {
                intercambiar(idx, mayor);
                idx = mayor;
            } else {
                break;
            }
        }
    }

    void insertar(const string& kmer, int conteo) {
        if (tamano >= capacidad) {
            int nuevaCap = capacidad * 2;
            NodoHeap* nuevo = new NodoHeap[nuevaCap];
            for (int i = 0; i < tamano; i++) {
                nuevo[i] = heap[i];
            }
            delete[] heap;
            heap = nuevo;
            capacidad = nuevaCap;
        }
        heap[tamano].kmer = kmer;
        heap[tamano].conteo = conteo;
        subir(tamano);
        tamano++;
    }

    // Llena el Max-Heap con todos los k-mers guardados en la Tabla Hash 2
    void construirDesdeTablaHash(const TablaHashDoble& tablaHash2) {
        tamano = 0;
        for (int i = 0; i < tablaHash2.capacidad; i++) {
            if (tablaHash2.tabla[i].ocupado) {
                insertar(tablaHash2.tabla[i].kmer, tablaHash2.tabla[i].conteo);
            }
        }
    }

    // Extrae y retorna la raiz (el k-mer de mayor frecuencia)
    NodoHeap extraerMaximo() {
        NodoHeap maximo = heap[0];
        heap[0] = heap[tamano - 1];
        tamano--;
        if (tamano > 0) {
            bajar(0);
        }
        return maximo;
    }
};
