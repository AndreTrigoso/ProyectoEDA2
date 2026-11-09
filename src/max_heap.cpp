#pragma once
#include <iostream>
#include <string>
#include "tablas_hash.cpp"

using namespace std;

// ============================================================================
// c) MAX-HEAP PARA RETORNAR LOS N PATRONES (K-MERS) MAS FRECUENTES
// ============================================================================

struct MaxHeapKmers {
    KmerConteo* arreglo;
    int tamano;
    int capacidad;

    MaxHeapKmers() {
        capacidad = 64;
        tamano = 0;
        arreglo = new KmerConteo[capacidad];
    }

    ~MaxHeapKmers() {
        delete[] arreglo;
    }

    // Compara si el elemento 'a' debe estar mas arriba que 'b' en el Max-Heap
    bool esMayor(const KmerConteo& a, const KmerConteo& b) const {
        if (a.conteo != b.conteo) {
            return a.conteo > b.conteo;
        }
        return a.kmer < b.kmer; // Desempate alfabetico
    }

    void intercambiar(int i, int j) {
        KmerConteo temp = arreglo[i];
        arreglo[i] = arreglo[j];
        arreglo[j] = temp;
    }

    // Sube un nodo si su frecuencia es mayor que la de su padre
    void subir(int idx) {
        while (idx > 0) {
            int padre = (idx - 1) / 2;
            if (esMayor(arreglo[idx], arreglo[padre])) {
                intercambiar(idx, padre);
                idx = padre;
            } else {
                break;
            }
        }
    }

    // Baja un nodo si alguno de sus hijos tiene mayor frecuencia
    void bajar(int idx) {
        while (true) {
            int hijoIzq = 2 * idx + 1;
            int hijoDer = 2 * idx + 2;
            int mayor = idx;

            if (hijoIzq < tamano && esMayor(arreglo[hijoIzq], arreglo[mayor])) {
                mayor = hijoIzq;
            }
            if (hijoDer < tamano && esMayor(arreglo[hijoDer], arreglo[mayor])) {
                mayor = hijoDer;
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
            KmerConteo* nuevoArr = new KmerConteo[nuevaCap];
            for (int i = 0; i < tamano; i++) {
                nuevoArr[i] = arreglo[i];
            }
            delete[] arreglo;
            arreglo = nuevoArr;
            capacidad = nuevaCap;
        }
        arreglo[tamano] = KmerConteo(kmer, conteo);
        subir(tamano);
        tamano++;
    }

    // Carga todos los k-mers ocupados desde la Tabla Hash 2 al Max-Heap
    void construirDesdeTablaHash(const TablaHashDoble& tablaHash2) {
        delete[] arreglo;
        capacidad = (tablaHash2.kmersUnicos > 0) ? tablaHash2.kmersUnicos : 16;
        arreglo = new KmerConteo[capacidad];
        tamano = 0;

        for (int i = 0; i < tablaHash2.capacidad; i++) {
            if (tablaHash2.tabla[i].ocupado) {
                insertar(tablaHash2.tabla[i].kmer, tablaHash2.tabla[i].conteo);
            }
        }
    }

    // Extrae la raiz (el k-mer con mayor frecuencia) y reordena el Max-Heap
    KmerConteo extraerMaximo() {
        if (tamano <= 0) {
            return KmerConteo("", 0);
        }
        KmerConteo maximo = arreglo[0];
        arreglo[0] = arreglo[tamano - 1];
        tamano--;
        if (tamano > 0) {
            bajar(0);
        }
        return maximo;
    }

    // Devuelve un arreglo con los N patrones mas frecuentes sin vaciar el heap original
    KmerConteo* obtenerTopN(int n, int& cantidadExtraida) const {
        if (n <= 0 || tamano == 0) {
            cantidadExtraida = 0;
            return nullptr;
        }

        cantidadExtraida = (n < tamano) ? n : tamano;

        // Copia auxiliar para extraer los N mayores
        MaxHeapKmers copia;
        delete[] copia.arreglo;
        copia.capacidad = tamano;
        copia.tamano = tamano;
        copia.arreglo = new KmerConteo[tamano];
        for (int i = 0; i < tamano; i++) {
            copia.arreglo[i] = arreglo[i];
        }

        KmerConteo* top = new KmerConteo[cantidadExtraida];
        for (int i = 0; i < cantidadExtraida; i++) {
            top[i] = copia.extraerMaximo();
        }
        return top;
    }
};
