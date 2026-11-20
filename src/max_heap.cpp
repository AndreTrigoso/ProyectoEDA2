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
    // funcion auxiliar para intercambiar dos nodos del arreglo del heap
    NodoHeap temp = heap[i];
    heap[i] = heap[j];
    heap[j] = temp;
}

void MaxHeapKmers::subir(int idx) {
    // heapify up: cuando metemos un k-mer al final lo subimos hasta que respete la propiedad de max-heap
    // mientras el indice sea mayor a 0 (o sea mientras no estemos ya en la raiz)
    while (idx > 0) {
        // en un arreglo 0-indexado el padre siempre esta en (idx - 1) / 2
        int padre = (idx - 1) / 2;
        comparacionesHeap++; // [GENERADO POR IA PARA CONTRASTACION]
        // si la frecuencia del hijo supera a la del padre, lo subimos intercambiando posiciones
        if (heap[idx].conteo > heap[padre].conteo) {
            intercambiar(idx, padre);
            idx = padre;
        } else {
            // si el padre ya es mayor o igual, el nodo quedo bien acomodado y paramos
            break;
        }
    }
}

void MaxHeapKmers::bajar(int idx) {
    // heapify down: cuando sacamos la raiz bajamos el nodo que quedo arriba hasta restaurar el orden
    while (true) {
        // en un arreglo 0-indexado los hijos de idx quedan en 2*idx + 1 y 2*idx + 2
        int izq = 2 * idx + 1;
        int der = 2 * idx + 2;
        int mayor = idx;

        // comparamos primero con el hijo izquierdo si existe dentro del tamano actual
        if (izq < tamano) {
            comparacionesHeap++; // [GENERADO POR IA PARA CONTRASTACION]
            if (heap[izq].conteo > heap[mayor].conteo) {
                mayor = izq;
            }
        }
        // luego comparamos con el hijo derecho para quedarnos con el mayor de los tres
        if (der < tamano) {
            comparacionesHeap++; // [GENERADO POR IA PARA CONTRASTACION]
            if (heap[der].conteo > heap[mayor].conteo) {
                mayor = der;
            }
        }

        // si alguno de los hijos era mas grande, intercambiamos y seguimos bajando desde ahi
        if (mayor != idx) {
            intercambiar(idx, mayor);
            idx = mayor;
        } else {
            // si el nodo ya es mayor o igual que sus dos hijos, el subarbol quedo valido
            break;
        }
    }
}

void MaxHeapKmers::insertar(const string& kmer, int conteo) {
    // si ya llenamos el arreglo dinamico, duplicamos su capacidad y copiamos los nodos actuales
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
    // colocamos el nuevo k-mer al final del arbol (primera posicion libre) y lo hacemos subir
    heap[tamano].kmer = kmer;
    heap[tamano].conteo = conteo;
    subir(tamano);
    tamano++;
}

void MaxHeapKmers::construirDesdeTablaHash(const TablaHashDoble& tablaHash2) {
    // reiniciamos el heap y recorremos los 4^k bins de la tabla hash 2
    tamano = 0;
    comparacionesHeap = 0; // [GENERADO POR IA PARA CONTRASTACION]
    // cada casilla ocupada trae un k-mer distinto con su frecuencia acumulada y se inserta en el heap
    for (int i = 0; i < tablaHash2.capacidad; i++) {
        if (tablaHash2.tabla[i].ocupado) {
            insertar(tablaHash2.tabla[i].kmer, tablaHash2.tabla[i].conteo);
        }
    }
}

NodoHeap MaxHeapKmers::extraerMaximo() {
    // proteccion por si alguien llama extraerMaximo con el heap vacio
    if (tamano <= 0) {
        return {"", 0};
    }
    // la frecuencia mas alta siempre esta en la raiz (heap[0])
    NodoHeap maximo = heap[0];
    // movemos el ultimo elemento del arreglo a la raiz para mantener el arbol completo
    heap[0] = heap[tamano - 1];
    tamano--;
    // si todavia quedan elementos, bajamos la nueva raiz para que el siguiente maximo suba a heap[0]
    if (tamano > 0) {
        bajar(0);
    }
    return maximo;
}
