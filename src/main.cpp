#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

// Inclusion directa de los archivos .cpp del proyecto (sin headers .hpp)
#include "parser.cpp"
#include "tablas_hash.cpp"
#include "max_heap.cpp"
#include "grafo_adn.cpp"

using namespace std;

// Funcion auxiliar para escapar texto al enviar JSON a la interfaz de Python
string escaparJson(const string& s) {
    ostringstream oss;
    for (int i = 0; i < (int)s.length(); i++) {
        char c = s[i];
        if (c == '"') oss << "\\\"";
        else if (c == '\\') oss << "\\\\";
        else if (c == '\n') oss << "\\n";
        else if (c == '\r') oss << "\\r";
        else if (c == '\t') oss << "\\t";
        else oss << c;
    }
    return oss.str();
}

void imprimirEstadoJson(const FastaParser& parser, const TablaHashLineal& ht1, bool ok, const string& msg) {
    ostringstream out;
    out << fixed << setprecision(6);
    out << "{";
    out << "\"ok\":" << (ok ? "true" : "false") << ",";
    out << "\"message\":\"" << escaparJson(msg) << "\",";
    out << "\"fastaCount\":" << parser.cantidad << ",";
    out << "\"ht1Capacity\":" << ht1.capacidad << ",";
    out << "\"ht1Count\":" << ht1.elementos << ",";
    out << "\"ht1LoadFactor\":" << ht1.obtenerLoadFactor() << ",";
    out << "\"ht1Collisions\":" << ht1.colisionesTotales << ",";
    out << "\"ht1Rehashes\":[";
    for (int i = 0; i < ht1.cantidadRehashes; i++) {
        if (i > 0) out << ",";
        out << "{\"oldCap\":" << ht1.historialRehash[i].capAnterior
            << ",\"newCap\":" << ht1.historialRehash[i].capNueva
            << ",\"elements\":" << ht1.historialRehash[i].elementos
            << ",\"loadFactor\":" << ht1.historialRehash[i].nuevoLoadFactor << "}";
    }
    out << "],";
    out << "\"sequences\":[";
    for (int i = 0; i < parser.cantidad; i++) {
        if (i > 0) out << ",";
        const DNASequence& seq = parser.lista[i];
        string preview = seq.sequence.length() > 36 ? (seq.sequence.substr(0, 36) + "...") : seq.sequence;
        out << "{\"index\":" << i
            << ",\"id\":\"" << escaparJson(seq.id) << "\""
            << ",\"length\":" << seq.sequence.length()
            << ",\"preview\":\"" << escaparJson(preview) << "\"}";
    }
    out << "]}";
    cout << out.str() << endl;
}

void imprimirGrafoJson(const GrafoADN& grafo, const string& etiqueta, const string& seqOriginal) {
    ostringstream out;
    out << fixed << setprecision(5);
    out << "{";
    out << "\"ok\":true,";
    out << "\"label\":\"" << escaparJson(etiqueta) << "\",";
    out << "\"seqLength\":" << seqOriginal.length() << ",";
    out << "\"beta\":" << grafo.beta << ",";
    out << "\"pairCount\":" << grafo.cantidadPares << ",";

    int limitePares = (grafo.cantidadPares < 300) ? grafo.cantidadPares : 300;
    out << "\"pairs\":[";
    for (int i = 0; i < limitePares; i++) {
        if (i > 0) out << ",";
        out << "{\"pos\":" << grafo.listaPares[i].posicion
            << ",\"pair\":\"" << grafo.listaPares[i].par << "\"}";
    }
    out << "],";

    int vActivos = grafo.cantidadVerticesActivos;
    out << "\"activeLabels\":[";
    for (int i = 0; i < vActivos; i++) {
        if (i > 0) out << ",";
        out << "\"" << grafo.verticesActivos[i] << "\"";
    }
    out << "],";

    out << "\"activeMatrix\":[";
    for (int i = 0; i < vActivos; i++) {
        if (i > 0) out << ",";
        out << "[";
        for (int j = 0; j < vActivos; j++) {
            if (j > 0) out << ",";
            out << grafo.matrizActiva[i][j];
        }
        out << "]";
    }
    out << "],";

    out << "\"canonicalLabels\":[";
    for (int i = 0; i < 16; i++) {
        if (i > 0) out << ",";
        out << "\"" << GrafoADN::indiceCanonicoAPar(i) << "\"";
    }
    out << "],";

    out << "\"canonicalMatrix\":[";
    for (int i = 0; i < 16; i++) {
        if (i > 0) out << ",";
        out << "[";
        for (int j = 0; j < 16; j++) {
            if (j > 0) out << ",";
            out << grafo.matrizCanonica[i][j];
        }
        out << "]";
    }
    out << "],";

    out << "\"w_TC_AT\":" << grafo.obtenerPeso("TC", "AT") << ",";
    out << "\"w_AT_CA\":" << grafo.obtenerPeso("AT", "CA");
    out << "}";
    cout << out.str() << endl;
}

// Modo puente para conectarse con la interfaz grafica en Python (gui/app.py)
int ejecutarModoServidor(const string& rutaFastaInicial) {
    FastaParser parser;
    TablaHashLineal ht1(10); // Tabla Hash 1 con Linear Probing (empieza en 10)
    TablaHashDoble ht2(3);   // Tabla Hash 2 con Double Hashing (4^n)
    MaxHeapKmers maxHeap;    // Max-Heap para los n patrones mas frecuentes
    GrafoADN grafo;          // Grafo de secuencias de ADN (Khan et al., 2024)

    if (!rutaFastaInicial.empty()) {
        if (parser.cargarArchivo(rutaFastaInicial)) {
            for (int i = 0; i < parser.cantidad; i++) {
                ht1.insertar(parser.lista[i]);
            }
        }
    }

    string linea;
    while (getline(cin, linea)) {
        if (linea.empty()) continue;
        istringstream iss(linea);
        string comando;
        iss >> comando;

        if (comando == "QUIT" || comando == "EXIT") {
            cout << "{\"ok\":true,\"message\":\"bye\"}" << endl;
            break;
        } else if (comando == "STATUS") {
            imprimirEstadoJson(parser, ht1, true, "Estado cargado");
        } else if (comando == "LOAD_FASTA") {
            string ruta;
            getline(iss >> ws, ruta);
            ht1.reiniciar();
            if (parser.cargarArchivo(ruta)) {
                for (int i = 0; i < parser.cantidad; i++) {
                    ht1.insertar(parser.lista[i]);
                }
                imprimirEstadoJson(parser, ht1, true, "Archivo FASTA cargado");
            } else {
                imprimirEstadoJson(parser, ht1, false, "Error al abrir el archivo FASTA");
            }
        } else if (comando == "SEARCH_SEQ") {
            string id;
            iss >> id;
            int sondeos = 0;
            const DNASequence* encontrado = ht1.buscar(id, &sondeos);
            ostringstream out;
            if (encontrado) {
                out << "{\"ok\":true,\"found\":true"
                    << ",\"id\":\"" << escaparJson(encontrado->id) << "\""
                    << ",\"length\":" << encontrado->sequence.length()
                    << ",\"comparisons\":" << sondeos
                    << ",\"sequence\":\"" << escaparJson(encontrado->sequence) << "\"}";
            } else {
                out << "{\"ok\":true,\"found\":false"
                    << ",\"id\":\"" << escaparJson(id) << "\""
                    << ",\"comparisons\":" << sondeos << "}";
            }
            cout << out.str() << endl;
        } else if (comando == "ANALYZE_KMERS") {
            string id;
            int k = 3;
            int topN = 10;
            iss >> id >> k >> topN;

            const DNASequence* seqObj = ht1.buscar(id);
            if (!seqObj) {
                cout << "{\"ok\":false,\"message\":\"Secuencia no encontrada en Tabla Hash 1\"}" << endl;
                continue;
            }

            ht2.configurarTamano(k);
            ht2.procesarSecuencia(seqObj->sequence);
            maxHeap.construirDesdeTablaHash(ht2);

            int extraidos = 0;
            KmerConteo* topArr = maxHeap.obtenerTopN(topN, extraidos);

            ostringstream out;
            out << fixed << setprecision(6);
            out << "{";
            out << "\"ok\":true,";
            out << "\"seqId\":\"" << escaparJson(seqObj->id) << "\",";
            out << "\"seqLength\":" << seqObj->sequence.length() << ",";
            out << "\"k\":" << ht2.longitudK << ",";
            out << "\"ht2Capacity\":" << ht2.capacidad << ",";
            out << "\"ht2Unique\":" << ht2.kmersUnicos << ",";
            out << "\"ht2TotalValid\":" << ht2.totalVentanasLeidas << ",";
            out << "\"ht2Ignored\":" << ht2.ventanasIgnoradas << ",";
            out << "\"ht2LoadFactor\":" << ht2.obtenerLoadFactor() << ",";
            out << "\"ht2Collisions\":" << ht2.colisionesTotales << ",";

            out << "\"topKmers\":[";
            for (int i = 0; i < extraidos; i++) {
                if (i > 0) out << ",";
                out << "{\"rank\":" << (i + 1)
                    << ",\"kmer\":\"" << escaparJson(topArr[i].kmer) << "\""
                    << ",\"count\":" << topArr[i].conteo << "}";
            }
            out << "],";
            delete[] topArr;

            out << "\"slots\":[";
            int mostrados = 0;
            for (int i = 0; i < ht2.capacidad && mostrados < 256; i++) {
                const CasillaKmer& s = ht2.tabla[i];
                if (ht2.capacidad <= 256 || s.ocupado) {
                    if (mostrados > 0) out << ",";
                    out << "{\"slot\":" << i
                        << ",\"occupied\":" << (s.ocupado ? "true" : "false")
                        << ",\"kmer\":\"" << (s.ocupado ? escaparJson(s.kmer) : "") << "\""
                        << ",\"count\":" << s.conteo
                        << ",\"probes\":" << s.sondeos << "}";
                    mostrados++;
                }
            }
            out << "]}";
            cout << out.str() << endl;
        } else if (comando == "SEARCH_KMER") {
            string kmer;
            iss >> kmer;
            for (int i = 0; i < (int)kmer.length(); i++) {
                if (kmer[i] >= 'a' && kmer[i] <= 'z') kmer[i] = kmer[i] - 'a' + 'A';
            }
            int sondeos = 0;
            int freq = ht2.buscar(kmer, &sondeos);
            ostringstream out;
            out << "{\"ok\":true"
                << ",\"kmer\":\"" << escaparJson(kmer) << "\""
                << ",\"found\":" << (freq > 0 ? "true" : "false")
                << ",\"count\":" << freq
                << ",\"comparisons\":" << sondeos
                << ",\"kLength\":" << ht2.longitudK
                << ",\"capacity\":" << ht2.capacidad << "}";
            cout << out.str() << endl;
        } else if (comando == "BUILD_GRAPH_ID") {
            string id;
            double beta = 1.0;
            iss >> id >> beta;
            const DNASequence* seqObj = ht1.buscar(id);
            if (!seqObj) {
                cout << "{\"ok\":false,\"message\":\"Secuencia no encontrada\"}" << endl;
                continue;
            }
            grafo.construirGrafo(seqObj->sequence, beta);
            imprimirGrafoJson(grafo, seqObj->id, seqObj->sequence);
        } else if (comando == "BUILD_GRAPH_RAW") {
            string seqPersonalizada;
            double beta = 1.0;
            iss >> seqPersonalizada >> beta;
            grafo.construirGrafo(seqPersonalizada, beta);
            imprimirGrafoJson(grafo, "Secuencia Personalizada", seqPersonalizada);
        } else {
            cout << "{\"ok\":false,\"message\":\"Comando desconocido\"}" << endl;
        }
    }
    return 0;
}

void imprimirMatrizConsola(const GrafoADN& grafo, bool mostrar16x16) {
    cout << fixed << setprecision(4);
    if (!mostrar16x16) {
        int n = grafo.cantidadVerticesActivos;
        cout << "\nMatriz de Adyacencia (" << n << "x" << n << ", beta = " << grafo.beta << "):\n      ";
        for (int j = 0; j < n; j++) {
            cout << setw(9) << grafo.verticesActivos[j];
        }
        cout << "\n";
        for (int i = 0; i < n; i++) {
            cout << setw(4) << grafo.verticesActivos[i] << "  ";
            for (int j = 0; j < n; j++) {
                cout << setw(9) << grafo.matrizActiva[i][j];
            }
            cout << "\n";
        }
    } else {
        cout << "\nMatriz de Adyacencia Canonica (16x16, beta = " << grafo.beta << "):\n     ";
        for (int j = 0; j < 16; j++) {
            cout << setw(8) << GrafoADN::indiceCanonicoAPar(j);
        }
        cout << "\n";
        for (int i = 0; i < 16; i++) {
            cout << setw(3) << GrafoADN::indiceCanonicoAPar(i) << "  ";
            for (int j = 0; j < 16; j++) {
                cout << setw(8) << grafo.matrizCanonica[i][j];
            }
            cout << "\n";
        }
    }
}

// e) Interfaz textual interactiva en consola C++
int ejecutarMenuConsola(const string& rutaFasta) {
    FastaParser parser;
    TablaHashLineal ht1(10);
    TablaHashDoble ht2(3);
    MaxHeapKmers maxHeap;
    GrafoADN grafo;

    cout << "==============================================================\n";
    cout << "         MINI SISTEMA DE ANALISIS GENOMICO (C++)              \n";
    cout << "==============================================================\n";

    if (parser.cargarArchivo(rutaFasta)) {
        for (int i = 0; i < parser.cantidad; i++) {
            ht1.insertar(parser.lista[i]);
        }
        cout << "[OK] Archivo cargado: " << rutaFasta << "\n";
        cout << "     Secuencias almacenadas: " << parser.cantidad << "\n";
        cout << "     Tabla Hash 1 (Linear Probing) -> Capacidad (primo m): "
             << ht1.capacidad << " | Load Factor: " << fixed << setprecision(4)
             << ht1.obtenerLoadFactor() << "\n";
    } else {
        cout << "[!] No se pudo abrir " << rutaFasta << "\n";
    }

    int opcion = -1;
    while (opcion != 0) {
        cout << "\n--- MENU PRINCIPAL ---\n";
        cout << "1. (a/b) Ver redimensionamientos de Tabla Hash 1 y buscar secuencia por ID\n";
        cout << "2. (b/c) Contar k-mers con Double Hashing (4^n) y ver Top-N con Max-Heap\n";
        cout << "3. (d)   Construir Grafo y Matriz de Adyacencia de una secuencia HIV\n";
        cout << "4. (d)   Probar ejemplo del PDF de Grafo (TCATTAGCAGTTAGCAGC, beta = 1)\n";
        cout << "5. (b)   Probar ejemplo corto del PDF de 3-mers (ATGGGATGTCTTG)\n";
        cout << "0. Salir\n";
        cout << "Seleccione una opcion: ";
        if (!(cin >> opcion)) break;

        if (opcion == 1) {
            cout << "\nHistorial de Redimensionamiento de Tabla Hash 1 (Inicio = 10):\n";
            for (int i = 0; i < ht1.cantidadRehashes; i++) {
                cout << "  Rehash #" << (i + 1) << ": " << ht1.historialRehash[i].capAnterior
                     << " -> Primo m = " << ht1.historialRehash[i].capNueva
                     << " (Elementos: " << ht1.historialRehash[i].elementos
                     << ", Nuevo LF: " << ht1.historialRehash[i].nuevoLoadFactor << ")\n";
            }
            cout << "\nIngrese ID de la secuencia (ej. MAC.US.x.239.M33262): ";
            string id;
            cin >> id;
            int sondeos = 0;
            const DNASequence* seq = ht1.buscar(id, &sondeos);
            if (seq) {
                cout << "  [Encontrado en " << sondeos << " sondeo(s)]\n";
                cout << "  ID: " << seq->id << " | Longitud: " << seq->sequence.length() << "\n";
                cout << "  Inicio de cadena: " << seq->sequence.substr(0, 60) << "...\n";
            } else {
                cout << "  [No encontrado tras " << sondeos << " sondeo(s)]\n";
            }
        } else if (opcion == 2) {
            string id;
            int k = 3, topN = 10;
            cout << "Ingrese ID de la secuencia (ej. MAC.US.x.239.M33262): ";
            cin >> id;
            const DNASequence* seq = ht1.buscar(id);
            if (!seq) {
                cout << "  [!] Secuencia no encontrada.\n";
                continue;
            }
            cout << "Ingrese longitud k del k-mer (ej. 3 para 3-mers, 6 para 6-mers): ";
            cin >> k;
            cout << "Ingrese cantidad n de patrones mas frecuentes a mostrar (Max-Heap): ";
            cin >> topN;

            ht2.configurarTamano(k);
            ht2.procesarSecuencia(seq->sequence);
            maxHeap.construirDesdeTablaHash(ht2);

            cout << "\n--- Tabla Hash 2 (Double Hashing) ---\n";
            cout << "  Tamano de la tabla (4^" << ht2.longitudK << "): " << ht2.capacidad << " bins\n";
            cout << "  K-mers distintos: " << ht2.kmersUnicos << "\n";
            cout << "  Load Factor:      " << ht2.obtenerLoadFactor() << "\n";
            cout << "  Colisiones:       " << ht2.colisionesTotales << "\n";

            int extraidos = 0;
            KmerConteo* top = maxHeap.obtenerTopN(topN, extraidos);
            cout << "\n--- Top " << extraidos << " k-mers mas frecuentes (Max-Heap) ---\n";
            for (int i = 0; i < extraidos; i++) {
                cout << "  #" << setw(2) << (i + 1) << ": " << top[i].kmer
                     << " -> " << top[i].conteo << " veces\n";
            }
            delete[] top;
        } else if (opcion == 3) {
            string id;
            double beta = 1.0;
            cout << "Ingrese ID de la secuencia (ej. MAC.US.x.239.M33262): ";
            cin >> id;
            const DNASequence* seq = ht1.buscar(id);
            if (!seq) {
                cout << "  [!] Secuencia no encontrada.\n";
                continue;
            }
            cout << "Ingrese valor de beta (ej. 1.0): ";
            cin >> beta;
            grafo.construirGrafo(seq->sequence, beta);
            imprimirMatrizConsola(grafo, true);
        } else if (opcion == 4) {
            string ejemplo = "TCATTAGCAGTTAGCAGC";
            grafo.construirGrafo(ejemplo, 1.0);
            cout << "\nCadena de ejemplo: " << ejemplo << " (beta = 1)\n";
            cout << "Posicion | 2-mer\n";
            for (int i = 0; i < grafo.cantidadPares; i++) {
                cout << setw(8) << grafo.listaPares[i].posicion << " | " << grafo.listaPares[i].par << "\n";
            }
            cout << fixed << setprecision(5);
            cout << "\nPesos del ejemplo del PDF:\n";
            cout << "  w(TC, AT) = " << grafo.obtenerPeso("TC", "AT") << "\n";
            cout << "  w(AT, CA) = " << grafo.obtenerPeso("AT", "CA") << "\n";
            imprimirMatrizConsola(grafo, false);
        } else if (opcion == 5) {
            string prefijo = "ATGGGATGTCTTG";
            ht2.configurarTamano(3);
            ht2.procesarSecuencia(prefijo);
            cout << "\nEjemplo del PDF (Prefijo de MAC.US.x.239.M33262): " << prefijo << "\n";
            cout << "Tamano Tabla Hash 2 (4^3): " << ht2.capacidad << " bins\n";
            cout << "Key (3-mers) | Valor (Frecuencia)\n";
            cout << "  ATG        | " << ht2.buscar("ATG") << "\n";
            cout << "  TGG        | " << ht2.buscar("TGG") << "\n";
            cout << "  GGG        | " << ht2.buscar("GGG") << "\n";
            cout << "  GGA        | " << ht2.buscar("GGA") << "\n";
            cout << "  GAT        | " << ht2.buscar("GAT") << "\n";
        }
    }
    return 0;
}

int main(int argc, char* argv[]) {
    string rutaFasta = "HIV2_ALL_2021_env_DNA.2000.fasta";
    bool modoServidor = false;

    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--server") {
            modoServidor = true;
            if (i + 1 < argc) {
                rutaFasta = argv[++i];
            }
        }
    }

    if (modoServidor) {
        return ejecutarModoServidor(rutaFasta);
    }
    return ejecutarMenuConsola(rutaFasta);
}
