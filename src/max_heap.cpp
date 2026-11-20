#include "../include/max_heap.hpp"

MaxHeapKmers::MaxHeapKmers() {
    capacidad = 64;
    tamano = 0;
    comparacionesHeap = 0; // [GENERADO POR IA PARA CONTRASTACION]
    heap = new NodoHeap[capacidad];
}

MaxHeapKmers::~MaxHeapKmers() {
    delete[] heap;
}

void MaxHeapKmers::intercambiar(int i, int j) {
    //funcion auxiliar para intercambio de nodos
    NodoHeap temp = heap[i];
    heap[i] = heap[j];
    heap[j] = temp;
}

void MaxHeapKmers::subir(int idx) {
    //heapofy yp
    // mientras que el indice del nodo sea mayor a 0
    while (idx > 0) {
        //obtenemos el padre
        int padre = (idx - 1) / 2;
        comparacionesHeap++; // [GENERADO POR IA PARA CONTRASTACION]
        if (heap[idx].conteo > heap[padre].conteo) {
            intercambiar(idx, padre);
            idx = padre;
        } else {
            break;
        }
    }
}

void MaxHeapKmers::bajar(int idx) {
    // heapify down
    while (true) {
        int izq = 2 * idx + 1;
        int der = 2 * idx + 2;
        int mayor = idx;

        if (izq < tamano) {
            comparacionesHeap++; // [GENERADO POR IA PARA CONTRASTACION]
            if (heap[izq].conteo > heap[mayor].conteo) {
                mayor = izq;
            }
        }
        if (der < tamano) {
            comparacionesHeap++; // [GENERADO POR IA PARA CONTRASTACION]
            if (heap[der].conteo > heap[mayor].conteo) {
                mayor = der;
            }
        }

        if (mayor != idx) {
            intercambiar(idx, mayor);
            idx = mayor;
        } else {
            break;
        }
    }
}

void MaxHeapKmers::insertar(const string& kmer, int conteo) {
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

void MaxHeapKmers::construirDesdeTablaHash(const TablaHashDoble& tablaHash2) {
    tamano = 0;
    comparacionesHeap = 0; // [GENERADO POR IA PARA CONTRASTACION]
    for (int i = 0; i < tablaHash2.capacidad; i++) {
        if (tablaHash2.tabla[i].ocupado) {
            insertar(tablaHash2.tabla[i].kmer, tablaHash2.tabla[i].conteo);
        }
    }
}

NodoHeap MaxHeapKmers::extraerMaximo() {
    NodoHeap maximo = heap[0];
    heap[0] = heap[tamano - 1];
    tamano--;
    if (tamano > 0) {
        bajar(0);
    }
    return maximo;
}
