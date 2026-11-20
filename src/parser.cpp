#include "../include/parser.hpp"
#include <fstream>
#include <sstream>
#include <cctype>

// ============================================================================
// INCISO a) PARSER DE SECUENCIAS ADN (FASTA y TXT) CON ARREGLO DINAMICO
// ============================================================================

bool esKmerValido(const string& kmer) {
    for (int i = 0; i < (int)kmer.length(); i++) {
        char c = kmer[i];
        if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
            return false;
        }
    }
    return true;
}

FastaParser::FastaParser() {
    capacidad = 20;
    cantidad = 0;
    lista = new DNASequence[capacidad];
}

FastaParser::~FastaParser() {
    delete[] lista;
}

void FastaParser::agrandarLista() {
    int nuevaCap = capacidad * 2;
    DNASequence* nuevaLista = new DNASequence[nuevaCap];
    for (int i = 0; i < cantidad; i++) {
        nuevaLista[i] = lista[i];
    }
    delete[] lista;
    lista = nuevaLista;
    capacidad = nuevaCap;
}

bool FastaParser::cargarArchivo(const string& rutaArchivo) {
    ifstream archivo(rutaArchivo.c_str());
    if (!archivo.is_open()) return false;

    cantidad = 0;
    string linea;

    DNASequence actual;
    actual.id = "";
    actual.sequence = "";

    while (getline(archivo, linea)) {
        // Limpiar saltos de linea de Windows ('\r') y espacios al final
        while (!linea.empty() && (linea[linea.length() - 1] == '\r' ||
                                  linea[linea.length() - 1] == '\n' ||
                                  linea[linea.length() - 1] == ' ')) {
            linea.erase(linea.length() - 1);
        }
        if (linea.empty()) continue;

        if (linea[0] == '>') {
            // Guardar la secuencia anterior si ya teniamos una acumulada
            if (actual.id != "" && actual.sequence != "") {
                if (cantidad >= capacidad) agrandarLista();
                lista[cantidad] = actual;
                cantidad++;
            }

            // Extraer el nuevo ID despues del simbolo '>'
            stringstream ss(linea.substr(1));
            ss >> actual.id;
            actual.sequence = "";
        }
        else if (actual.id != "") {
            // Formato FASTA: acumular solo nucleotidos validos en mayuscula
            for (int i = 0; i < (int)linea.length(); i++) {
                char c = toupper((unsigned char)linea[i]);
                if (c == 'A' || c == 'T' || c == 'C' || c == 'G') {
                    actual.sequence += c;
                }
            }
        }
        else {
            // Formato alternativo .txt en una sola linea: "ID SECUENCIA"
            stringstream ss(linea);
            DNASequence dna;
            string parte;
            ss >> dna.id;
            dna.sequence = "";

            while (ss >> parte) {
                for (int i = 0; i < (int)parte.length(); i++) {
                    char c = toupper((unsigned char)parte[i]);
                    if (c == 'A' || c == 'T' || c == 'C' || c == 'G') {
                        dna.sequence += c;
                    }
                }
            }

            if (dna.id != "" && dna.sequence != "") {
                if (cantidad >= capacidad) agrandarLista();
                lista[cantidad] = dna;
                cantidad++;
            }
        }
    }

    // Guardar la ultima secuencia del archivo FASTA
    if (actual.id != "" && actual.sequence != "") {
        if (cantidad >= capacidad) agrandarLista();
        lista[cantidad] = actual;
        cantidad++;
    }

    archivo.close();
    return cantidad > 0;
}