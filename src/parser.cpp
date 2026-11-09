#pragma once
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

// a) Struct sugerido en el documento para almacenar cada secuencia con su ID
struct DNASequence {
    string id;       // Ej: MAC.US.x.239.M33262
    string sequence; // Ej: ATGGGATGTCTTGGGAA...
};

// Verifica si un caracter es uno de los 4 nucleotidos principales (A, T, C, G)
bool esNucleotidoValido(char c) {
    return (c == 'A' || c == 'T' || c == 'C' || c == 'G');
}

// Verifica si un k-mer completo solo contiene A, T, C, G
bool esKmerValido(const string& kmer) {
    if (kmer.empty()) return false;
    for (int i = 0; i < (int)kmer.length(); i++) {
        if (!esNucleotidoValido(kmer[i])) return false;
    }
    return true;
}

// Parser para leer el archivo FASTA de HIV
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

    void redimensionarLista() {
        int nuevaCap = capacidad * 2;
        DNASequence* nuevaLista = new DNASequence[nuevaCap];
        for (int i = 0; i < cantidad; i++) {
            nuevaLista[i] = lista[i];
        }
        delete[] lista;
        lista = nuevaLista;
        capacidad = nuevaCap;
    }

    void limpiar() {
        cantidad = 0;
    }

    bool cargarArchivo(const string& rutaArchivo) {
        ifstream archivo(rutaArchivo.c_str());
        if (!archivo.is_open()) {
            return false;
        }

        limpiar();
        string linea;
        string idActual = "";
        string seqActual = "";
        bool leyendoSecuencia = false;

        while (getline(archivo, linea)) {
            // Quitar saltos de linea o espacios al final
            while (!linea.empty() && (linea.back() == '\r' || linea.back() == '\n' || linea.back() == ' ')) {
                linea.pop_back();
            }
            if (linea.empty()) continue;

            if (linea[0] == '>') {
                // Guardar la secuencia anterior antes de empezar la nueva
                if (leyendoSecuencia && !idActual.empty()) {
                    if (cantidad >= capacidad) {
                        redimensionarLista();
                    }
                    lista[cantidad].id = idActual;
                    lista[cantidad].sequence = seqActual;
                    cantidad++;
                }
                idActual = linea.substr(1); // Quitar el caracter '>'
                seqActual = "";
                leyendoSecuencia = true;
            } else if (leyendoSecuencia) {
                for (int i = 0; i < (int)linea.length(); i++) {
                    char c = linea[i];
                    if (c >= 'a' && c <= 'z') c = c - 'a' + 'A';
                    if (c != ' ' && c != '\t') {
                        seqActual += c;
                    }
                }
            }
        }

        // Guardar la ultima secuencia del archivo
        if (leyendoSecuencia && !idActual.empty()) {
            if (cantidad >= capacidad) {
                redimensionarLista();
            }
            lista[cantidad].id = idActual;
            lista[cantidad].sequence = seqActual;
            cantidad++;
        }

        archivo.close();
        return cantidad > 0;
    }
};
