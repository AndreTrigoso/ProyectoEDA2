#include "../include/parser.hpp"
#include <fstream>
#include <sstream>
#include <cctype>

// ============================================================================
// INCISO a) PARSER DE SECUENCIAS ADN (FASTA y TXT) CON ARREGLO DINAMICO
// ============================================================================

bool esKmerValido(const string& kmer) {
    if (kmer.empty()) return false;
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
    bool primeraLinea = true;

    DNASequence actual;
    actual.id = "";
    actual.sequence = "";

    while (getline(archivo, linea)) {
        if (primeraLinea) {
            primeraLinea = false;
            if (linea.size() >= 3 &&
                static_cast<unsigned char>(linea[0]) == 0xEF &&
                static_cast<unsigned char>(linea[1]) == 0xBB &&
                static_cast<unsigned char>(linea[2]) == 0xBF) {
                linea.erase(0, 3);
            }
        }

        // Limpiar saltos de linea de Windows ('\r') y espacios al final
        while (!linea.empty() && (linea[linea.length() - 1] == '\r' ||
                                  linea[linea.length() - 1] == '\n' ||
                                  linea[linea.length() - 1] == ' ' ||
                                  linea[linea.length() - 1] == '\t')) {
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
            // Formato FASTA: conservar todas las letras en mayuscula (incluidas
            // bases ambiguas IUPAC como N, R, Y) para no corromper la fase de lectura
            for (int i = 0; i < (int)linea.length(); i++) {
                unsigned char uc = static_cast<unsigned char>(linea[i]);
                if (isalpha(uc)) {
                    char c = static_cast<char>(toupper(uc));
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
                    unsigned char uc = static_cast<unsigned char>(parte[i]);
                    if (isalpha(uc)) {
                        char c = static_cast<char>(toupper(uc));
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