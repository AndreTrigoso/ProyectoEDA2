#pragma once
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

// a) Struct sugerido en el documento para cada secuencia
struct DNASequence {
    string id;       // MAC.US.x.239.M33262
    string sequence; // ACTG...
};

// Verifica que el k-mer solo tenga nucleotidos validos (A, T, C, G)
bool esKmerValido(const string& kmer) {
    for (int i = 0; i < (int)kmer.length(); i++) {
        char c = kmer[i];
        if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
            return false;
        }
    }
    return true;
}

// Parser para leer las secuencias del archivo FASTA
struct FastaParser {
    DNASequence* lista;
    int cantidad;
    int capacidad;

    FastaParser() {
        capacidad = 20;
        cantidad = 0;
        lista = new DNASequence[capacidad];
    }

    ~FastaParser() {
        delete[] lista;
    }

    void agrandarLista() {
        int nuevaCap = capacidad * 2;
        DNASequence* nuevaLista = new DNASequence[nuevaCap];
        for (int i = 0; i < cantidad; i++) {
            nuevaLista[i] = lista[i];
        }
        delete[] lista;
        lista = nuevaLista;
        capacidad = nuevaCap;
    }

    bool cargarArchivo(const string& rutaArchivo) {
        ifstream archivo(rutaArchivo.c_str());
        if (!archivo.is_open()) return false;

        cantidad = 0;
        string linea;
        string idActual = "";
        string seqActual = "";

        while (getline(archivo, linea)) {
            while (!linea.empty() && (linea.back() == '\r' || linea.back() == '\n' || linea.back() == ' ')) {
                linea.pop_back();
            }
            if (linea.empty()) continue;

            if (linea[0] == '>') {
                if (!idActual.empty()) {
                    if (cantidad >= capacidad) agrandarLista();
                    lista[cantidad].id = idActual;
                    lista[cantidad].sequence = seqActual;
                    cantidad++;
                }
                idActual = linea.substr(1); // Quitar '>'
                seqActual = "";
            } else {
                seqActual += linea;
            }
        }

        // Guardar la ultima secuencia del archivo
        if (!idActual.empty()) {
            if (cantidad >= capacidad) agrandarLista();
            lista[cantidad].id = idActual;
            lista[cantidad].sequence = seqActual;
            cantidad++;
        }

        archivo.close();
        return cantidad > 0;
    }
};
