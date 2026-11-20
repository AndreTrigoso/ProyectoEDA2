#include "../include/tablas_hash.hpp"

// ============================================================================
// FUNCIONES MATEMATICAS BASICAS
// ============================================================================

bool esPrimo(int n) {
    if (n <= 1) return false;
    for (int i = 2; i * i <= n; i++) {
        if (n % i == 0) return false;
    }
    return true;
}

int siguientePrimo(int numero) {
    while (!esPrimo(numero)) {
        numero++;
    }
    return numero;
}

int potenciaDeCuatro(int n) {
    int resultado = 1;
    for (int i = 0; i < n; i++) {
        resultado *= 4;
    }
    return resultado;
}

// ============================================================================
// b.1) TABLA HASH 1: LINEAR PROBING (ID -> Secuencia)
// ============================================================================

CasillaSeq::CasillaSeq() {
    ocupado = false;
}

TablaHashLineal::TablaHashLineal(int tamInicial) {
    capacidad = tamInicial;
    elementos = 0;
    tabla = new CasillaSeq[capacidad];

    // [GENERADO POR IA PARA CONTRASTACION EN LA DISCUSION]
    totalSondas = 0;
    colisiones = 0;
    maxSondas = 0;
    numRehashes = 0;
    ultimasSondasBusqueda = 0;
    enRehash = false;
}

TablaHashLineal::~TablaHashLineal() {
    delete[] tabla;
}

unsigned long TablaHashLineal::funcionHash(const string& id) const {
    unsigned long hash = 0;
    for (int i = 0; i < (int)id.length(); i++) {
        hash = (hash * 31) + id[i];
    }
    return hash % capacidad;
}

double TablaHashLineal::obtenerLoadFactor() const {
    return (double)elementos / (double)capacidad;
}

void TablaHashLineal::redimensionar(int nuevoTamanoPrimo) {
    CasillaSeq* tablaVieja = tabla;
    int capVieja = capacidad;

    capacidad = nuevoTamanoPrimo;
    elementos = 0;
    tabla = new CasillaSeq[capacidad];

    // [GENERADO POR IA PARA CONTRASTACION] Evita doble conteo durante reubicacion interna
    numRehashes++;
    enRehash = true;
    for (int i = 0; i < capVieja; i++) {
        if (tablaVieja[i].ocupado) {
            insertar(tablaVieja[i].dato);
        }
    }
    enRehash = false;
    delete[] tablaVieja;
}

void TablaHashLineal::insertar(const DNASequence& seq) {
    // Si al insertar el load factor supera 1, redimensionar al siguiente
    // primo m que baje el load factor a 0.5 (m = siguientePrimo(2 * nuevoTotal))
    if ((double)(elementos + 1) / capacidad > 1.0) {
        int primoM = siguientePrimo((elementos + 1) * 2);
        redimensionar(primoM);
    }

    unsigned long posInicial = funcionHash(seq.id);

    // Linear Probing: avanza de 1 en 1 si la casilla esta ocupada
    for (int i = 0; i < capacidad; i++) {
        int idx = (posInicial + i) % capacidad;

        if (!tabla[idx].ocupado) {
            // [GENERADO POR IA PARA CONTRASTACION EN LA DISCUSION]
            if (!enRehash) {
                int sondasActuales = i + 1;
                totalSondas += sondasActuales;
                colisiones += i;
                if (sondasActuales > maxSondas) maxSondas = sondasActuales;
            }
            tabla[idx].dato = seq;
            tabla[idx].ocupado = true;
            elementos++;
            return;
        }
        if (tabla[idx].dato.id == seq.id) {
            if (!enRehash) {
                int sondasActuales = i + 1;
                totalSondas += sondasActuales;
                colisiones += i;
                if (sondasActuales > maxSondas) maxSondas = sondasActuales;
            }
            tabla[idx].dato.sequence = seq.sequence;
            return;
        }
    }
}

const DNASequence* TablaHashLineal::buscar(const string& id) const {
    unsigned long posInicial = funcionHash(id);
    ultimasSondasBusqueda = 0; // [GENERADO POR IA PARA CONTRASTACION]

    for (int i = 0; i < capacidad; i++) {
        int idx = (posInicial + i) % capacidad;
        ultimasSondasBusqueda = i + 1; // [GENERADO POR IA PARA CONTRASTACION]

        if (!tabla[idx].ocupado) {
            return nullptr; // no existe
        }
        if (tabla[idx].dato.id == id) {
            return &(tabla[idx].dato);
        }
    }
    return nullptr;
}

void TablaHashLineal::reiniciar() {
    delete[] tabla;
    capacidad = 10;
    elementos = 0;
    tabla = new CasillaSeq[capacidad];

    // [GENERADO POR IA PARA CONTRASTACION EN LA DISCUSION]
    totalSondas = 0;
    colisiones = 0;
    maxSondas = 0;
    numRehashes = 0;
    ultimasSondasBusqueda = 0;
    enRehash = false;
}

// ============================================================================
// b.2) TABLA HASH 2: DOUBLE HASHING (k-mer -> Frecuencia, tamano 4^n)
// ============================================================================

CasillaKmer::CasillaKmer() {
    kmer = "";
    conteo = 0;
    ocupado = false;
    sondas = 0; // [GENERADO POR IA PARA CONTRASTACION]
}

TablaHashDoble::TablaHashDoble(int k) {
    tabla = nullptr;
    configurarTamano(k);
}

TablaHashDoble::~TablaHashDoble() {
    delete[] tabla;
}

void TablaHashDoble::configurarTamano(int nuevoK) {
    longitudK = nuevoK;
    capacidad = potenciaDeCuatro(longitudK);
    kmersUnicos = 0;

    // [GENERADO POR IA PARA CONTRASTACION EN LA DISCUSION]
    totalSondas = 0;
    colisiones = 0;
    maxSondas = 0;
    sondasLinearProbing = 0;
    colisionesLinearProbing = 0;
    maxSondasLinearProbing = 0;
    ultimasSondasBusqueda = 0;

    delete[] tabla;
    tabla = new CasillaKmer[capacidad];
}

unsigned long TablaHashDoble::hash1(const string& kmer) const {
    unsigned long h = 0;
    for (int i = 0; i < (int)kmer.length(); i++) {
        h = (h * 37) + kmer[i];
    }
    return h % capacidad;
}

unsigned long TablaHashDoble::hash2(const string& kmer) const {
    unsigned long h = 0;
    for (int i = 0; i < (int)kmer.length(); i++) {
        h = (h * 53) + kmer[i];
    }
    return ((h % (capacidad / 2)) * 2) + 1;
}

void TablaHashDoble::insertarKmer(const string& kmer) {
    if ((int)kmer.length() != longitudK || !esKmerValido(kmer)) return;

    unsigned long h1 = hash1(kmer);
    unsigned long h2 = hash2(kmer);

    // Double Hashing: (h1 + i * h2) % capacidad
    for (int i = 0; i < capacidad; i++) {
        int idx = (h1 + i * h2) % capacidad;

        if (!tabla[idx].ocupado) {
            // [GENERADO POR IA PARA CONTRASTACION EN LA DISCUSION]
            int sondasActuales = i + 1;
            totalSondas += sondasActuales;
            colisiones += i;
            if (sondasActuales > maxSondas) maxSondas = sondasActuales;

            tabla[idx].kmer = kmer;
            tabla[idx].conteo = 1;
            tabla[idx].ocupado = true;
            tabla[idx].sondas = sondasActuales;
            kmersUnicos++;
            return;
        }
        if (tabla[idx].kmer == kmer) {
            // [GENERADO POR IA PARA CONTRASTACION EN LA DISCUSION]
            int sondasActuales = i + 1;
            totalSondas += sondasActuales;
            colisiones += i;
            if (sondasActuales > maxSondas) maxSondas = sondasActuales;

            tabla[idx].conteo++; // Si ya existe, suma 1 a su frecuencia
            return;
        }
    }
}

void TablaHashDoble::procesarSecuencia(const string& secuencia) {
    configurarTamano(longitudK);
    int largo = (int)secuencia.length();

    for (int pos = 0; pos <= largo - longitudK; pos++) {
        string subcadena = secuencia.substr(pos, longitudK);
        insertarKmer(subcadena);
    }

    // [GENERADO POR IA PARA CONTRASTACION EN LA DISCUSION]
    contrastarConLinearProbing(secuencia);
}

// ----------------------------------------------------------------------------
// [GENERADO POR IA PARA CONTRASTACION EN LA DISCUSION]
// Prompt utilizado: "Implementa una simulacion auxiliar de Linear Probing
// usando exactamente la misma funcion hash1 y la misma capacidad 4^k para
// medir sondasLinearProbing, colisionesLinearProbing y maxSondasLinearProbing
// y poder contrastar Double Hashing vs Linear Probing en la discusion."
// ----------------------------------------------------------------------------
void TablaHashDoble::contrastarConLinearProbing(const string& secuencia) {
    sondasLinearProbing = 0;
    colisionesLinearProbing = 0;
    maxSondasLinearProbing = 0;

    CasillaKmer* tablaLP = new CasillaKmer[capacidad];
    int largo = (int)secuencia.length();

    for (int pos = 0; pos <= largo - longitudK; pos++) {
        string kmer = secuencia.substr(pos, longitudK);
        if ((int)kmer.length() != longitudK || !esKmerValido(kmer)) continue;

        unsigned long h1 = hash1(kmer);
        for (int i = 0; i < capacidad; i++) {
            int idx = (h1 + i) % capacidad;
            if (!tablaLP[idx].ocupado) {
                int sondasActuales = i + 1;
                sondasLinearProbing += sondasActuales;
                colisionesLinearProbing += i;
                if (sondasActuales > maxSondasLinearProbing) {
                    maxSondasLinearProbing = sondasActuales;
                }
                tablaLP[idx].kmer = kmer;
                tablaLP[idx].conteo = 1;
                tablaLP[idx].ocupado = true;
                break;
            }
            if (tablaLP[idx].kmer == kmer) {
                int sondasActuales = i + 1;
                sondasLinearProbing += sondasActuales;
                colisionesLinearProbing += i;
                if (sondasActuales > maxSondasLinearProbing) {
                    maxSondasLinearProbing = sondasActuales;
                }
                tablaLP[idx].conteo++;
                break;
            }
        }
    }
    delete[] tablaLP;
}

int TablaHashDoble::buscar(const string& kmer) const {
    ultimasSondasBusqueda = 0; // [GENERADO POR IA PARA CONTRASTACION]
    if ((int)kmer.length() != longitudK || !esKmerValido(kmer)) return 0;

    unsigned long h1 = hash1(kmer);
    unsigned long h2 = hash2(kmer);

    for (int i = 0; i < capacidad; i++) {
        int idx = (h1 + i * h2) % capacidad;
        ultimasSondasBusqueda = i + 1; // [GENERADO POR IA PARA CONTRASTACION]

        if (!tabla[idx].ocupado) return 0;
        if (tabla[idx].kmer == kmer) return tabla[idx].conteo;
    }
    return 0;
}
