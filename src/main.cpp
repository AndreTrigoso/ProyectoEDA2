#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <limits>
#include <cmath>
#include <cctype>

#include "../include/parser.hpp"
#include "../include/tablas_hash.hpp"
#include "../include/max_heap.hpp"
#include "../include/grafo_adn.hpp"

using namespace std;

// Lee un entero desde cin validando todo el token:
// retorna 1 si es valido, 0 si tiene formato invalido, -1 si hay EOF
int leerEntero(int& valor) {
    string token;
    if (!(cin >> token)) {
        return -1;
    }
    istringstream iss(token);
    char extra;
    if (!(iss >> valor) || (iss >> extra)) {
        if (cin.peek() != '\n' && cin.peek() != EOF) {
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
        }
        return 0;
    }
    return 1;
}

// Lee un double desde cin validando todo el token:
// retorna 1 si es valido, 0 si tiene formato invalido, -1 si hay EOF
int leerDouble(double& valor) {
    string token;
    if (!(cin >> token)) {
        return -1;
    }
    istringstream iss(token);
    char extra;
    if (!(iss >> valor) || (iss >> extra)) {
        if (cin.peek() != '\n' && cin.peek() != EOF) {
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
        }
        return 0;
    }
    return 1;
}

// Imprime la matriz de adyacencia del grafo en consola
void imprimirMatriz(const GrafoADN& grafo) {
    int n = grafo.numVertices;
    if (n <= 0) {
        cout << "\nEl grafo no tiene vertices validos.\n";
        return;
    }
    cout << fixed << setprecision(4);
    cout << "\nMatriz de Adyacencia (" << n << "x" << n << ", beta = " << grafo.beta << "):\n      ";
    for (int j = 0; j < n; j++) {
        cout << " " << setw(9) << grafo.vertices[j];
    }
    cout << "\n";
    for (int i = 0; i < n; i++) {
        cout << setw(4) << grafo.vertices[i] << "  ";
        for (int j = 0; j < n; j++) {
            cout << " " << setw(9) << grafo.matriz[i][j];
        }
        cout << "\n";
    }
}

// Modo servidor para conectar con la interfaz grafica en Python (gui/app.py)
int ejecutarModoServidor(const string& rutaFasta) {
    FastaParser parser;
    TablaHashLineal ht1(10);
    TablaHashDoble ht2(3);
    MaxHeapKmers maxHeap;
    GrafoADN grafo;

    if (parser.cargarArchivo(rutaFasta)) {
        for (int i = 0; i < parser.cantidad; i++) {
            ht1.insertar(parser.lista[i]);
        }
    }

    string linea;
    while (getline(cin, linea)) {
        if (linea.empty()) continue;
        istringstream iss(linea);
        string cmd;
        iss >> cmd;

        if (cmd == "QUIT") {
            cout << "{\"ok\":true}" << endl;
            break;
        } else if (cmd == "STATUS" || cmd == "LOAD_FASTA") {
            if (cmd == "LOAD_FASTA") {
                string nuevaRuta;
                getline(iss >> ws, nuevaRuta);
                ht1.reiniciar();
                if (parser.cargarArchivo(nuevaRuta)) {
                    for (int i = 0; i < parser.cantidad; i++) {
                        ht1.insertar(parser.lista[i]);
                    }
                }
            }
            ostringstream out;
            out << fixed << setprecision(4);
            // [GENERADO POR IA PARA CONTRASTACION] Usa ht1.colisiones reales
            out << "{\"ok\":true,\"fastaCount\":" << parser.cantidad
                << ",\"ht1Capacity\":" << ht1.capacidad
                << ",\"ht1Count\":" << ht1.elementos
                << ",\"ht1LoadFactor\":" << ht1.obtenerLoadFactor()
                << ",\"ht1Collisions\":" << ht1.colisiones
                << ",\"ht1Rehashes\":[],\"sequences\":[";
            for (int i = 0; i < parser.cantidad; i++) {
                if (i > 0) out << ",";
                string prev = parser.lista[i].sequence.substr(0, 36) + "...";
                out << "{\"index\":" << i << ",\"id\":\"" << parser.lista[i].id
                    << "\",\"length\":" << parser.lista[i].sequence.length()
                    << ",\"preview\":\"" << prev << "\"}";
            }
            out << "]}";
            cout << out.str() << endl;
        } else if (cmd == "SEARCH_SEQ") {
            string id;
            iss >> id;
            const DNASequence* seq = ht1.buscar(id);
            if (seq) {
                // [GENERADO POR IA PARA CONTRASTACION] Expone ultimasSondasBusqueda reales
                cout << "{\"ok\":true,\"found\":true,\"id\":\"" << seq->id
                     << "\",\"length\":" << seq->sequence.length()
                     << ",\"comparisons\":" << ht1.ultimasSondasBusqueda
                     << ",\"sequence\":\"" << seq->sequence << "\"}" << endl;
            } else {
                cout << "{\"ok\":true,\"found\":false,\"id\":\"" << id
                     << "\",\"comparisons\":" << ht1.ultimasSondasBusqueda << "}" << endl;
            }
        } else if (cmd == "ANALYZE_KMERS") {
            string id;
            int k = 3, topN = 10;
            iss >> id >> k >> topN;
            if (k < 1 || k > MAX_K_VALIDO) {
                cout << "{\"ok\":false,\"message\":\"k fuera de rango (1..10)\"}" << endl;
                continue;
            }
            const DNASequence* seq = ht1.buscar(id);
            if (!seq) {
                cout << "{\"ok\":false,\"message\":\"No encontrado\"}" << endl;
                continue;
            }
            ht2.configurarTamano(k);
            ht2.procesarSecuencia(seq->sequence);
            maxHeap.construirDesdeTablaHash(ht2);

            int totalValidos = ht2.totalSondas - ht2.colisiones;
            ostringstream out;
            out << fixed << setprecision(4);
            // [GENERADO POR IA PARA CONTRASTACION] Expone ht2.colisiones y metricas reales
            out << "{\"ok\":true,\"seqId\":\"" << seq->id
                << "\",\"seqLength\":" << seq->sequence.length()
                << ",\"k\":" << ht2.longitudK
                << ",\"ht2Capacity\":" << ht2.capacidad
                << ",\"ht2Unique\":" << ht2.kmersUnicos
                << ",\"ht2TotalValid\":" << totalValidos
                << ",\"ht2LoadFactor\":" << (double)ht2.kmersUnicos / ht2.capacidad
                << ",\"ht2Collisions\":" << ht2.colisiones << ",\"topKmers\":[";

            int limite = (topN < maxHeap.tamano) ? topN : maxHeap.tamano;
            for (int i = 0; i < limite; i++) {
                if (i > 0) out << ",";
                NodoHeap maximo = maxHeap.extraerMaximo();
                out << "{\"rank\":" << (i + 1) << ",\"kmer\":\"" << maximo.kmer
                    << "\",\"count\":" << maximo.conteo << "}";
            }
            out << "],\"slots\":[";
            int mostrados = 0;
            for (int i = 0; i < ht2.capacidad && mostrados < 256; i++) {
                if (ht2.capacidad <= 256 || ht2.tabla[i].ocupado) {
                    if (mostrados > 0) out << ",";
                    out << "{\"slot\":" << i
                        << ",\"occupied\":" << (ht2.tabla[i].ocupado ? "true" : "false")
                        << ",\"kmer\":\"" << ht2.tabla[i].kmer << "\""
                        << ",\"count\":" << ht2.tabla[i].conteo
                        << ",\"probes\":" << ht2.tabla[i].sondas << "}";
                    mostrados++;
                }
            }
            out << "]}";
            cout << out.str() << endl;
        } else if (cmd == "SEARCH_KMER") {
            string kmer;
            iss >> kmer;
            int freq = ht2.buscar(kmer);
            // [GENERADO POR IA PARA CONTRASTACION] Expone ht2.ultimasSondasBusqueda reales
            cout << "{\"ok\":true,\"kmer\":\"" << kmer
                 << "\",\"found\":" << (freq > 0 ? "true" : "false")
                 << ",\"count\":" << freq << ",\"comparisons\":" << ht2.ultimasSondasBusqueda
                 << ",\"kLength\":" << ht2.longitudK
                 << ",\"capacity\":" << ht2.capacidad << "}" << endl;
        } else if (cmd == "BUILD_GRAPH_ID" || cmd == "BUILD_GRAPH_RAW") {
            string entrada;
            double beta = 1.0;
            iss >> entrada >> beta;
            if (!(beta > 0.0) || !isfinite(beta)) {
                beta = 1.0;
            }
            string secuencia = entrada;
            string etiqueta = "Secuencia Personalizada";
            if (cmd == "BUILD_GRAPH_ID") {
                const DNASequence* seq = ht1.buscar(entrada);
                if (!seq) {
                    cout << "{\"ok\":false}" << endl;
                    continue;
                }
                secuencia = seq->sequence;
                etiqueta = seq->id;
            }
            grafo.construirGrafo(secuencia, beta);

            ostringstream out;
            out << fixed << setprecision(5);
            out << "{\"ok\":true,\"label\":\"" << etiqueta
                << "\",\"seqLength\":" << secuencia.length()
                << ",\"beta\":" << grafo.beta
                << ",\"pairCount\":" << grafo.cantidadPares << ",\"pairs\":[";
            int maxP = (grafo.cantidadPares < 300) ? grafo.cantidadPares : 300;
            for (int i = 0; i < maxP; i++) {
                if (i > 0) out << ",";
                out << "{\"pos\":" << grafo.listaPares[i].posicion
                    << ",\"pair\":\"" << grafo.listaPares[i].par << "\"}";
            }
            out << "],\"activeLabels\":[";
            for (int i = 0; i < grafo.numVertices; i++) {
                if (i > 0) out << ",";
                out << "\"" << grafo.vertices[i] << "\"";
            }
            out << "],\"activeMatrix\":[";
            for (int i = 0; i < grafo.numVertices; i++) {
                if (i > 0) out << ",";
                out << "[";
                for (int j = 0; j < grafo.numVertices; j++) {
                    if (j > 0) out << ",";
                    out << grafo.matriz[i][j];
                }
                out << "]";
            }
            out << "],\"canonicalLabels\":[],\"canonicalMatrix\":[],"
                << "\"w_TC_AT\":" << grafo.obtenerPeso("TC", "AT") << ","
                << "\"w_AT_CA\":" << grafo.obtenerPeso("AT", "CA") << "}";
            cout << out.str() << endl;
        }
    }
    return 0;
}

// e) Menu interactivo por consola en C++
int main(int argc, char* argv[]) {
    string rutaFasta = "HIV2_ALL_2021_env_DNA.2000.fasta";
    if (argc > 1 && string(argv[1]) == "--server") {
        if (argc > 2) rutaFasta = argv[2];
        return ejecutarModoServidor(rutaFasta);
    }

    FastaParser parser;
    TablaHashLineal ht1(10); // Empieza con tamano 10
    TablaHashDoble ht2(3);
    MaxHeapKmers maxHeap;
    GrafoADN grafo;
    bool ht2Cargada = false;
    string seqActualHt2 = "";

    if (parser.cargarArchivo(rutaFasta)) {
        for (int i = 0; i < parser.cantidad; i++) {
            ht1.insertar(parser.lista[i]);
        }
        cout << "Secuencias cargadas: " << parser.cantidad << "\n";
        cout << "Tamano final Tabla Hash 1 (primo m): " << ht1.capacidad
             << " | Load Factor: " << ht1.obtenerLoadFactor() << "\n";
        // [GENERADO POR IA PARA CONTRASTACION EN LA DISCUSION]
        // Prompt utilizado: "Imprime las metricas de sondas, colisiones, maximo
        // de sondas y rehashes de TablaHashLineal y la contrastacion Double
        // Hashing vs Linear Probing y MaxHeap en consola."
        cout << fixed << setprecision(4);
        cout << "[Contrastacion HT1 - Linear Probing] Rehashes: " << ht1.numRehashes
             << " | Sondas totales: " << ht1.totalSondas
             << " | Promedio: " << (double)ht1.totalSondas / parser.cantidad
             << " | Colisiones: " << ht1.colisiones
             << " | Max sondas: " << ht1.maxSondas << "\n";
    } else {
        cout << "[ERROR] No se pudo abrir o parsear el archivo FASTA: " << rutaFasta << "\n";
        cout << "        Las opciones 4 y 5 (ejemplos del PDF) siguen disponibles.\n";
    }

    int opcion = -1;
    while (opcion != 0) {
        cout << "\n--- MENU ---\n";
        cout << "1. Buscar secuencia por ID en Tabla Hash 1 (Linear Probing)\n";
        cout << "2. Contar k-mers en Tabla Hash 2 (Double Hashing 4^n) y Top-N con Max-Heap\n";
        cout << "3. Construir Grafo y Matriz de Adyacencia de una secuencia HIV\n";
        cout << "4. Probar ejemplo del PDF de Grafo (TCATTAGCAGTTAGCAGC, beta = 1)\n";
        cout << "5. Probar ejemplo del PDF de 3-mers (ATGGGATGTCTTG)\n";
        cout << "6. Buscar un k-mer especifico en Tabla Hash 2 (Double Hashing)\n";
        cout << "0. Salir\n";
        cout << "Opcion: ";

        int estadoOpt = leerEntero(opcion);
        if (estadoOpt == -1) break;
        if (estadoOpt == 0) {
            cout << "Entrada invalida. Ingrese un numero entre 0 y 6.\n";
            opcion = -1;
            continue;
        }

        if (opcion == 0) {
            break;
        } else if (opcion == 1) {
            string id;
            cout << "ID a buscar (ej. MAC.US.x.239.M33262): ";
            if (!(cin >> id)) break;
            const DNASequence* seq = ht1.buscar(id);
            if (seq) {
                cout << "Encontrado! Longitud: " << seq->sequence.length() << "\n";
                cout << "Prefijo: " << seq->sequence.substr(0, 50) << "...\n";
                // [GENERADO POR IA PARA CONTRASTACION EN LA DISCUSION]
                cout << "[Contrastacion] Sondas en busqueda HT1: " << ht1.ultimasSondasBusqueda << "\n";
            } else {
                cout << "No encontrado (sondas: " << ht1.ultimasSondasBusqueda << ").\n";
            }
        } else if (opcion == 2) {
            string id;
            int k = 0, n = 0;
            cout << "ID de secuencia (ej. MAC.US.x.239.M33262): ";
            if (!(cin >> id)) break;
            const DNASequence* seq = ht1.buscar(id);
            if (!seq) {
                cout << "Secuencia no encontrada.\n";
                continue;
            }
            cout << "Longitud de k-mer (1-" << MAX_K_VALIDO << ", ej. 3 o 6): ";
            int estadoK = leerEntero(k);
            if (estadoK == -1) break;
            if (estadoK == 0 || k < 1 || k > MAX_K_VALIDO) {
                cout << "Longitud k invalida. Debe ser un entero entre 1 y " << MAX_K_VALIDO << ".\n";
                continue;
            }
            cout << "Cantidad n de patrones mas frecuentes (Max-Heap): ";
            int estadoN = leerEntero(n);
            if (estadoN == -1) break;
            if (estadoN == 0 || n < 1) {
                cout << "Cantidad n invalida. Debe ser un entero mayor o igual a 1.\n";
                continue;
            }

            ht2.configurarTamano(k);
            ht2.procesarSecuencia(seq->sequence);
            maxHeap.construirDesdeTablaHash(ht2);
            ht2Cargada = true;
            seqActualHt2 = seq->id;

            cout << "\nTamano Tabla Hash 2 (4^" << k << "): " << ht2.capacidad << " bins\n";
            cout << "K-mers distintos guardados: " << ht2.kmersUnicos << "\n";
            cout << "\nTop " << n << " k-mers mas frecuentes:\n";
            for (int i = 0; i < n && maxHeap.tamano > 0; i++) {
                NodoHeap top = maxHeap.extraerMaximo();
                cout << "  #" << (i + 1) << ": " << top.kmer << " -> " << top.conteo << "\n";
            }
            // [GENERADO POR IA PARA CONTRASTACION EN LA DISCUSION]
            int ventanasValidas = ht2.totalSondas - ht2.colisiones;
            cout << fixed << setprecision(4);
            cout << "\n[Contrastacion HT2 (M = " << ht2.capacidad << ", ventanas = " << ventanasValidas << ")]\n";
            cout << "  Double Hashing -> Sondas totales: " << ht2.totalSondas
                 << " | Promedio: " << (ventanasValidas > 0 ? (double)ht2.totalSondas / ventanasValidas : 0.0)
                 << " | Colisiones: " << ht2.colisiones
                 << " | Max sondas: " << ht2.maxSondas << "\n";
            cout << "  Linear Probing -> Sondas totales: " << ht2.sondasLinearProbing
                 << " | Promedio: " << (ventanasValidas > 0 ? (double)ht2.sondasLinearProbing / ventanasValidas : 0.0)
                 << " | Colisiones: " << ht2.colisionesLinearProbing
                 << " | Max sondas: " << ht2.maxSondasLinearProbing << "\n";
            cout << "[Contrastacion Max-Heap] Comparaciones realizadas (build + Top-" << n << "): "
                 << maxHeap.comparacionesHeap << "\n";
        } else if (opcion == 3) {
            string id;
            double beta = 0.0;
            cout << "ID de secuencia (ej. MAC.US.x.239.M33262): ";
            if (!(cin >> id)) break;
            const DNASequence* seq = ht1.buscar(id);
            if (!seq) {
                cout << "Secuencia no encontrada.\n";
                continue;
            }
            cout << "Valor de beta (> 0, ej. 1.0): ";
            int estadoBeta = leerDouble(beta);
            if (estadoBeta == -1) break;
            if (estadoBeta == 0 || !(beta > 0.0) || !isfinite(beta)) {
                cout << "Valor de beta invalido. Debe ser un numero positivo mayor que 0.\n";
                continue;
            }
            grafo.construirGrafo(seq->sequence, beta);
            imprimirMatriz(grafo);
        } else if (opcion == 4) {
            grafo.construirGrafo("TCATTAGCAGTTAGCAGC", 1.0);
            cout << "\nPosicion | 2-mer\n";
            for (int i = 0; i < grafo.cantidadPares; i++) {
                cout << setw(8) << grafo.listaPares[i].posicion << " | " << grafo.listaPares[i].par << "\n";
            }
            cout << fixed << setprecision(5);
            cout << "\nw(TC, AT) = " << grafo.obtenerPeso("TC", "AT") << "\n";
            cout << "w(AT, CA) = " << grafo.obtenerPeso("AT", "CA") << "\n";
            imprimirMatriz(grafo);
        } else if (opcion == 5) {
            ht2.configurarTamano(3);
            ht2.procesarSecuencia("ATGGGATGTCTTG");
            ht2Cargada = true;
            seqActualHt2 = "Ejemplo PDF (ATGGGATGTCTTG)";
            cout << "\nCadena: ATGGGATGTCTTG (Tabla Hash 2 = " << ht2.capacidad << " bins)\n";
            cout << "ATG: " << ht2.buscar("ATG") << "\n";
            cout << "TGG: " << ht2.buscar("TGG") << "\n";
            cout << "GGG: " << ht2.buscar("GGG") << "\n";
            cout << "GGA: " << ht2.buscar("GGA") << "\n";
            cout << "GAT: " << ht2.buscar("GAT") << "\n";
        } else if (opcion == 6) {
            if (!ht2Cargada) {
                cout << "Tabla Hash 2 aun no tiene datos. Ejecute primero la opcion 2 o la opcion 5.\n";
                continue;
            }
            cout << "Tabla Hash 2 activa -> Secuencia: " << seqActualHt2
                 << " | k = " << ht2.longitudK
                 << " | Capacidad = " << ht2.capacidad << " bins\n";
            string kmer;
            cout << "K-mer a buscar (longitud " << ht2.longitudK << ", solo A/T/C/G): ";
            if (!(cin >> kmer)) break;

            string kmerMayus = "";
            for (int i = 0; i < (int)kmer.length(); i++) {
                kmerMayus += static_cast<char>(toupper(static_cast<unsigned char>(kmer[i])));
            }
            if ((int)kmerMayus.length() != ht2.longitudK) {
                cout << "Longitud invalida: el k-mer debe tener exactamente " << ht2.longitudK << " nucleotidos.\n";
                continue;
            }
            if (!esKmerValido(kmerMayus)) {
                cout << "K-mer invalido: solo se permiten los nucleotidos A, T, C y G.\n";
                continue;
            }
            int freq = ht2.buscar(kmer);
            if (freq > 0) {
                cout << "Encontrado! Frecuencia de " << kmerMayus << ": " << freq << "\n";
            } else {
                cout << "No encontrado (frecuencia = 0, ausente o posible nulomero en esta secuencia).\n";
            }
            // [GENERADO POR IA PARA CONTRASTACION EN LA DISCUSION]
            cout << "[Contrastacion] Sondas en busqueda HT2: " << ht2.ultimasSondasBusqueda << "\n";
        } else {
            cout << "Opcion invalida. Ingrese un numero entre 0 y 6.\n";
        }
    }
    return 0;
}
