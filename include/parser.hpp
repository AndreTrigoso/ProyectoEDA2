#ifndef PARSER_HPP
#define PARSER_HPP

#include <string>

using namespace std;

// a) Struct sugerido en el documento para cada secuencia
struct DNASequence {
    string id;       // MAC.US.x.239.M33262
    string sequence; // ACTG...
};

// Verifica que el k-mer solo tenga nucleotidos validos (A, T, C, G)
bool esKmerValido(const string& kmer);

// Parser para leer las secuencias del archivo FASTA
struct FastaParser {
    DNASequence* lista;
    int cantidad;
    int capacidad;

    FastaParser();
    ~FastaParser();

    void agrandarLista();
    bool cargarArchivo(const string& rutaArchivo);
};

#endif
