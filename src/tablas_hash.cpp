#include "../include/tablas_hash.hpp"
#include <cctype>

// ============================================================================
// FUNCIONES MATEMATICAS BASICAS
// ============================================================================

bool esPrimo(int n) {
    // los numeros <= 1 no son primos por definicion
    if (n <= 1) return false;
    // solo hace falta probar divisores hasta sqrt(n); usamos i <= n / i para evitar overflow en i*i
    for (int i = 2; i <= n / i; i++) {
        // si divide exacto entonces no es primo
        if (n % i == 0) return false;
    }
    return true;
}

int siguientePrimo(int numero) {
    // vamos avanzando de 1 en 1 desde el numero base hasta dar con el primer primo
    while (!esPrimo(numero)) {
        numero++;
    }
    return numero;
}

int potenciaDeCuatro(int n) {
    // caso base: 4^0 = 1
    if (n <= 0) return 1;
    // si n > 15 se desborda un int de 32 bits con signo, asi que devolvemos 0 como aviso
    if (n > 15) return 0;
    // multiplicamos por 4 n veces para calcular el numero exacto de bins (4^n)
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
    // si nos pasan un tamano invalido arrancamos en 10 como pide el enunciado
    capacidad = (tamInicial > 0) ? tamInicial : 10;
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
    // hash polinomial simple: multiplicamos por 31 y sumamos cada byte del id
    // usamos unsigned long porque su aritmetica modular es bien definida en c++
    unsigned long hash = 0;
    for (int i = 0; i < (int)id.length(); i++) {
        hash = (hash * 31UL) + static_cast<unsigned char>(id[i]);
    }
    // ajustamos el valor para que caiga dentro del arreglo [0, capacidad - 1]
    return hash % static_cast<unsigned long>(capacidad);
}

double TablaHashLineal::obtenerLoadFactor() const {
    return (double)elementos / (double)capacidad;
}

void TablaHashLineal::redimensionar(int nuevoTamanoPrimo) {
    // guardamos el puntero y tamano viejos para poder mover todas las secuencias a la tabla nueva
    CasillaSeq* tablaVieja = tabla;
    int capVieja = capacidad;

    // creamos la nueva tabla mas grande (de tamano primo m) y reiniciamos el conteo
    capacidad = nuevoTamanoPrimo;
    elementos = 0;
    tabla = new CasillaSeq[capacidad];

    // [GENERADO POR IA PARA CONTRASTACION] Evita doble conteo durante reubicacion interna
    numRehashes++;
    enRehash = true;
    // recorremos toda la tabla vieja y reinsertamos solo las casillas que estaban ocupadas
    // es necesario reinsertarlas porque al cambiar capacidad cambia su indice modulo capacidad
    for (int i = 0; i < capVieja; i++) {
        if (tablaVieja[i].ocupado) {
            insertar(tablaVieja[i].dato);
        }
    }
    enRehash = false;
    // liberamos la memoria del arreglo viejo para no tener fugas
    delete[] tablaVieja;
}

void TablaHashLineal::insertar(const DNASequence& seq) {
    // si meter un elemento nuevo pasaria el load factor de 1.0, primero revisamos si el id ya existia
    // asi evitamos redimensionar por gusto cuando solo estamos actualizando una secuencia existente
    if ((double)(elementos + 1) / capacidad > 1.0) {
        unsigned long posCheck = funcionHash(seq.id);
        for (int i = 0; i < capacidad; i++) {
            int idx = (int)((posCheck + (unsigned long)i) % (unsigned long)capacidad);
            if (!tabla[idx].ocupado) {
                // encontramos un hueco libre, o sea el id es nuevo y si toca hacer rehash
                break;
            }
            if (tabla[idx].dato.id == seq.id) {
                // el id ya estaba en la tabla, solo actualizamos su secuencia y salimos sin rehash
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

        // como el id es nuevo y el load factor pasaria de 1.0, buscamos el siguiente primo >= 2*(elementos+1)
        // ese primo garantiza que el nuevo load factor caiga a <= 0.5
        int primoM = siguientePrimo((elementos + 1) * 2);
        redimensionar(primoM);
    }

    // calculamos la casilla ideal donde deberia caer este id
    unsigned long posInicial = funcionHash(seq.id);

    // linear probing: si la casilla ideal esta ocupada avanzamos de 1 en 1 (dando la vuelta con %)
    for (int i = 0; i < capacidad; i++) {
        int idx = (int)((posInicial + (unsigned long)i) % (unsigned long)capacidad);

        // si hallamos una casilla libre, guardamos la secuencia ahi y aumentamos elementos
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
        // si nos topamos con el mismo id, actualizamos su secuencia sin duplicar la clave
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
    // arrancamos en la misma casilla inicial que usaria insertar
    unsigned long posInicial = funcionHash(id);
    ultimasSondasBusqueda = 0; // [GENERADO POR IA PARA CONTRASTACION]

    // vamos avanzando linealmente de 1 en 1 mientras encontremos casillas ocupadas
    for (int i = 0; i < capacidad; i++) {
        int idx = (int)((posInicial + (unsigned long)i) % (unsigned long)capacidad);
        ultimasSondasBusqueda = i + 1; // [GENERADO POR IA PARA CONTRASTACION]

        // si topamos una casilla vacia es seguro que el id nunca se inserto (no borramos claves aqui)
        if (!tabla[idx].ocupado) {
            return nullptr; // no existe
        }
        // si el id coincide devolvemos el puntero directo al registro encontrado
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
    capacidad = 0;
    longitudK = 0;
    kmersUnicos = 0;
    configurarTamano(k);
}

TablaHashDoble::~TablaHashDoble() {
    delete[] tabla;
}

void TablaHashDoble::configurarTamano(int nuevoK) {
    // acotamos k entre 1 y MAX_K_VALIDO para evitar division por cero o quedarnos sin RAM
    if (nuevoK < 1) nuevoK = 1;
    if (nuevoK > MAX_K_VALIDO) nuevoK = MAX_K_VALIDO;

    longitudK = nuevoK;
    int nuevaCap = potenciaDeCuatro(longitudK);
    kmersUnicos = 0;

    // [GENERADO POR IA PARA CONTRASTACION EN LA DISCUSION]
    totalSondas = 0;
    colisiones = 0;
    maxSondas = 0;
    sondasLinearProbing = 0;
    colisionesLinearProbing = 0;
    maxSondasLinearProbing = 0;
    ultimasSondasBusqueda = 0;

    // si la tabla ya tiene justo el tamano 4^k pedido, solo limpiamos las casillas sin pedir memoria otra vez
    if (tabla != nullptr && capacidad == nuevaCap) {
        for (int i = 0; i < capacidad; i++) {
            tabla[i] = CasillaKmer();
        }
    } else {
        // si cambio k, liberamos el arreglo anterior y reservamos exactamente 4^k bins
        delete[] tabla;
        capacidad = nuevaCap;
        tabla = new CasillaKmer[capacidad];
    }
}

unsigned long TablaHashDoble::hash1(const string& kmer) const {
    // primer hash polinomial con base prima 37 para ubicar el bin inicial
    unsigned long long h = 0;
    for (int i = 0; i < (int)kmer.length(); i++) {
        h = (h * 37ULL) + (unsigned char)kmer[i];
    }
    return (unsigned long)(h % (unsigned long long)capacidad);
}

unsigned long TablaHashDoble::hash2(const string& kmer) const {
    // segundo hash con primo 53 para calcular el tamano del salto cuando hay colision
    unsigned long long h = 0;
    for (int i = 0; i < (int)kmer.length(); i++) {
        h = (h * 53ULL) + (unsigned char)kmer[i];
    }
    // detalle importante: como capacidad = 4^k es potencia de 2, cualquier numero impar es coprimo con capacidad
    // multiplicar por 2 y sumar 1 garantiza que el paso siempre sea impar y >= 1,
    // asi en k=3 (64 bins) o cualquier k recorremos todos los bins sin caer en ciclos cortos
    return (unsigned long)(((h % (unsigned long long)(capacidad / 2)) * 2ULL) + 1ULL);
}

void TablaHashDoble::insertarKmer(const string& kmer) {
    // descartamos subcadenas que no midan k o que tengan letras ambiguas (ej. N, R, Y)
    if ((int)kmer.length() != longitudK || !esKmerValido(kmer)) return;

    unsigned long h1 = hash1(kmer);
    unsigned long h2 = hash2(kmer);

    // double hashing: en vez de movernos de 1 en 1, damos saltos de tamano h2 -> (h1 + i*h2) % capacidad
    // usamos unsigned long long para que i * h2 nunca haga overflow
    for (int i = 0; i < capacidad; i++) {
        int idx = (int)(((unsigned long long)h1 + (unsigned long long)i * (unsigned long long)h2) % (unsigned long long)capacidad);

        // si el bin esta vacio, es la primera vez que vemos este k-mer: lo registramos con conteo = 1
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
        // si ya estaba guardado ese mismo k-mer en este bin, solo le sumamos 1 a su frecuencia
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
    // reiniciamos la tabla con el tamano 4^longitudK antes de contar la secuencia
    configurarTamano(longitudK);
    int largo = (int)secuencia.length();

    // ventana deslizante con paso de 1 caracter: saca cada subcadena de largo k desde pos hasta largo - k
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
            int idx = (int)(((unsigned long long)h1 + (unsigned long long)i) % (unsigned long long)capacidad);
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

    // pasamos la consulta a mayusculas para que buscar "atg" funcione igual que "ATG"
    string kmerMayus = "";
    for (int i = 0; i < (int)kmer.length(); i++) {
        kmerMayus += static_cast<char>(toupper(static_cast<unsigned char>(kmer[i])));
    }
    // si no mide longitudK o tiene caracteres que no son A/T/C/G, no puede estar guardado
    if ((int)kmerMayus.length() != longitudK || !esKmerValido(kmerMayus)) return 0;

    // repetimos exactamente la misma secuencia de saltos (h1 + i*h2) % capacidad que uso insertarKmer
    unsigned long h1 = hash1(kmerMayus);
    unsigned long h2 = hash2(kmerMayus);

    for (int i = 0; i < capacidad; i++) {
        int idx = (int)(((unsigned long long)h1 + (unsigned long long)i * (unsigned long long)h2) % (unsigned long long)capacidad);
        ultimasSondasBusqueda = i + 1; // [GENERADO POR IA PARA CONTRASTACION]

        // si encontramos una casilla vacia antes de hallar el k-mer, significa que nunca aparecio (conteo 0)
        if (!tabla[idx].ocupado) return 0;
        // si coincide la clave, devolvemos cuantas veces aparecio en la secuencia
        if (tabla[idx].kmer == kmerMayus) return tabla[idx].conteo;
    }
    return 0;
}
