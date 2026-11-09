import json
import math
import os
import subprocess
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, ttk


class GenomicEngineClient:
    """Cliente que se comunica con el motor en C++ (genomic_engine.exe --server)."""

    def __init__(self, project_root: str):
        self.project_root = project_root
        self.exe_path = os.path.join(project_root, "genomic_engine.exe")
        self.fasta_path = os.path.join(project_root, "HIV2_ALL_2021_env_DNA.2000.fasta")
        self.proc = None
        self._start_engine()

    def _start_engine(self):
        if not os.path.exists(self.exe_path):
            build_script = os.path.join(self.project_root, "build.bat")
            if os.path.exists(build_script):
                subprocess.run([build_script], cwd=self.project_root, shell=True, check=False)

        if not os.path.exists(self.exe_path):
            raise FileNotFoundError(
                f"No se encontro el ejecutable C++ en:\n{self.exe_path}\n"
                "Ejecute build.bat primero para compilar el motor."
            )

        args = [self.exe_path, "--server"]
        if os.path.exists(self.fasta_path):
            args.append(self.fasta_path)

        creationflags = subprocess.CREATE_NO_WINDOW if sys.platform == "win32" else 0
        self.proc = subprocess.Popen(
            args,
            cwd=self.project_root,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            bufsize=1,
            creationflags=creationflags,
        )

    def send_command(self, command_line: str) -> dict:
        if self.proc is None or self.proc.poll() is not None:
            self._start_engine()
        self.proc.stdin.write(command_line.strip() + "\n")
        self.proc.stdin.flush()
        response_line = self.proc.stdout.readline()
        if not response_line:
            return {"ok": False, "message": "Sin respuesta del motor C++"}
        try:
            return json.loads(response_line)
        except json.JSONDecodeError as exc:
            return {"ok": False, "message": f"Error decodificando JSON: {exc}"}

    def close(self):
        if self.proc and self.proc.poll() is None:
            try:
                self.proc.stdin.write("QUIT\n")
                self.proc.stdin.flush()
                self.proc.terminate()
            except Exception:
                pass


class GenomicAnalyzerApp(tk.Tk):
    """Interfaz Grafica en Python para el Mini Sistema de Analisis Genomico (Punto e)."""

    def __init__(self, project_root: str):
        super().__init__()
        self.project_root = project_root
        self.title("Mini Sistema de Análisis Genómico — HIV-2 (Motor C++ + GUI Python)")
        self.geometry("1280x820")
        self.minsize(1100, 720)
        self.configure(bg="#0f172a")

        try:
            self.engine = GenomicEngineClient(project_root)
        except Exception as exc:
            messagebox.showerror("Error al iniciar motor C++", str(exc))
            self.destroy()
            return

        self.sequences_meta = []
        self.current_status = {}
        self.current_graph_data = None

        self._setup_styles()
        self._build_header()
        self._build_notebook()
        self._refresh_initial_data()

        self.protocol("WM_DELETE_WINDOW", self._on_close)

    def _setup_styles(self):
        style = ttk.Style(self)
        style.theme_use("clam")

        bg_dark = "#0f172a"
        bg_card = "#1e293b"
        fg_main = "#f8fafc"
        accent = "#38bdf8"
        accent_active = "#0284c7"

        style.configure(".", background=bg_dark, foreground=fg_main, font=("Segoe UI", 10))
        style.configure("TFrame", background=bg_dark)
        style.configure("Card.TFrame", background=bg_card, relief="flat")
        style.configure("TLabel", background=bg_dark, foreground=fg_main)
        style.configure("Card.TLabel", background=bg_card, foreground=fg_main)
        style.configure("Title.TLabel", font=("Segoe UI", 15, "bold"), foreground="#38bdf8")
        style.configure("Subtitle.TLabel", font=("Segoe UI", 9), foreground="#94a3b8")
        style.configure("MetricVal.TLabel", background=bg_card, foreground="#38bdf8", font=("Segoe UI", 15, "bold"))
        style.configure("MetricLbl.TLabel", background=bg_card, foreground="#cbd5e1", font=("Segoe UI", 9))

        style.configure(
            "TNotebook",
            background=bg_dark,
            borderwidth=0,
        )
        style.configure(
            "TNotebook.Tab",
            background="#1e293b",
            foreground="#cbd5e1",
            padding=(16, 8),
            font=("Segoe UI", 10, "bold"),
        )
        style.map(
            "TNotebook.Tab",
            background=[("selected", "#0284c7"), ("active", "#334155")],
            foreground=[("selected", "#ffffff"), ("active", "#ffffff")],
        )

        style.configure(
            "Accent.TButton",
            background=accent_active,
            foreground="#ffffff",
            padding=(12, 6),
            font=("Segoe UI", 10, "bold"),
            borderwidth=0,
        )
        style.map(
            "Accent.TButton",
            background=[("active", accent), ("pressed", "#0369a1")],
        )

        style.configure(
            "Treeview",
            background="#1e293b",
            foreground="#f8fafc",
            fieldbackground="#1e293b",
            rowheight=26,
            font=("Consolas", 10),
        )
        style.configure(
            "Treeview.Heading",
            background="#334155",
            foreground="#38bdf8",
            font=("Segoe UI", 9, "bold"),
        )
        style.map("Treeview", background=[("selected", "#0284c7")])

    def _build_header(self):
        header = ttk.Frame(self, padding=(16, 12))
        header.pack(fill="x")

        left = ttk.Frame(header)
        left.pack(side="left", fill="x", expand=True)
        ttk.Label(
            left,
            text="MINI SISTEMA DE ANÁLISIS GENÓMICO (C++ CORE)",
            style="Title.TLabel",
        ).pack(anchor="w")
        ttk.Label(
            left,
            text="a) Parser FASTA  |  b) Hash Table 1 (Linear Probing) & Hash Table 2 (Double Hashing 4^n)  |  c) Max-Heap Top-N  |  d) Grafo de ADN (Khan et al., 2024)",
            style="Subtitle.TLabel",
        ).pack(anchor="w")

        right = ttk.Frame(header)
        right.pack(side="right")
        ttk.Button(
            right,
            text="Cargar Otro Archivo FASTA...",
            style="Accent.TButton",
            command=self._on_load_fasta_dialog,
        ).pack(side="right")

    def _build_notebook(self):
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill="both", expand=True, padx=14, pady=(0, 14))

        self.tab_ht1 = ttk.Frame(self.notebook, padding=12)
        self.tab_ht2_heap = ttk.Frame(self.notebook, padding=12)
        self.tab_graph = ttk.Frame(self.notebook, padding=12)

        self.notebook.add(
            self.tab_ht1,
            text="  1. (a/b) Parser FASTA & Tabla Hash 1 (Linear Probing)  ",
        )
        self.notebook.add(
            self.tab_ht2_heap,
            text="  2. (b/c) Tabla Hash 2 (Double Hashing 4^n) & Max-Heap  ",
        )
        self.notebook.add(
            self.tab_graph,
            text="  3. (d) Grafo de Secuencias de ADN (Khan et al., 2024)  ",
        )

        self._build_tab_ht1()
        self._build_tab_ht2_heap()
        self._build_tab_graph()

    def _create_metric_card(self, parent, title: str, initial_val: str = "0"):
        card = ttk.Frame(parent, style="Card.TFrame", padding=(14, 10))
        val_lbl = ttk.Label(card, text=initial_val, style="MetricVal.TLabel")
        val_lbl.pack(anchor="w")
        ttk.Label(card, text=title, style="MetricLbl.TLabel").pack(anchor="w")
        return card, val_lbl

    # =========================================================================
    # PESTAÑA 1: PARSER FASTA & TABLA HASH 1 (LINEAR PROBING)
    # =========================================================================
    def _build_tab_ht1(self):
        metrics_row = ttk.Frame(self.tab_ht1)
        metrics_row.pack(fill="x", pady=(0, 10))

        c1, self.lbl_ht1_count = self._create_metric_card(metrics_row, "Secuencias en FASTA")
        c2, self.lbl_ht1_cap = self._create_metric_card(metrics_row, "Tamaño Tabla Hash 1 (Primo m)")
        c3, self.lbl_ht1_lf = self._create_metric_card(metrics_row, "Load Factor Actual")
        c4, self.lbl_ht1_rehashes = self._create_metric_card(metrics_row, "Redimensionamientos (10 → m)")
        c5, self.lbl_ht1_collisions = self._create_metric_card(metrics_row, "Colisiones Totales (Linear Probing)")

        for card in (c1, c2, c3, c4, c5):
            card.pack(side="left", fill="x", expand=True, padx=4)

        paned = ttk.PanedWindow(self.tab_ht1, orient="horizontal")
        paned.pack(fill="both", expand=True)

        # Panel izquierdo: Lista de secuencias y buscador por ID
        left_frame = ttk.Frame(paned, padding=(0, 0, 8, 0))
        paned.add(left_frame, weight=3)

        search_bar = ttk.Frame(left_frame)
        search_bar.pack(fill="x", pady=(0, 8))
        ttk.Label(search_bar, text="Buscar ID en Tabla Hash 1:").pack(side="left", padx=(0, 6))
        self.entry_search_id = ttk.Entry(search_bar, font=("Consolas", 10), width=28)
        self.entry_search_id.pack(side="left", padx=(0, 6))
        self.entry_search_id.insert(0, "MAC.US.x.239.M33262")
        self.entry_search_id.bind("<Return>", lambda _e: self._on_search_ht1())
        ttk.Button(
            search_bar,
            text="Buscar (Linear Probing)",
            style="Accent.TButton",
            command=self._on_search_ht1,
        ).pack(side="left")

        self.lbl_search_result = ttk.Label(
            left_frame,
            text="Seleccione o busque un identificador para inspeccionar su secuencia.",
            foreground="#38bdf8",
        )
        self.lbl_search_result.pack(anchor="w", pady=(0, 6))

        cols = ("idx", "id", "length", "preview")
        self.tree_sequences = ttk.Treeview(left_frame, columns=cols, show="headings", height=12)
        self.tree_sequences.heading("idx", text="#")
        self.tree_sequences.heading("id", text="ID Secuencia (Key)")
        self.tree_sequences.heading("length", text="Longitud (pb)")
        self.tree_sequences.heading("preview", text="Prefijo Cadena ADN (Valor)")
        self.tree_sequences.column("idx", width=45, anchor="center")
        self.tree_sequences.column("id", width=210, anchor="w")
        self.tree_sequences.column("length", width=95, anchor="center")
        self.tree_sequences.column("preview", width=260, anchor="w")
        self.tree_sequences.pack(fill="both", expand=True)
        self.tree_sequences.bind("<<TreeviewSelect>>", self._on_select_sequence_tree)

        # Panel derecho: Historial de Rehashing + Detalle de secuencia
        right_frame = ttk.Frame(paned, padding=(8, 0, 0, 0))
        paned.add(right_frame, weight=2)

        ttk.Label(
            right_frame,
            text="Historial de Redimensionamiento (Tamaño Inicial = 10 → Siguiente Primo m con LF ≈ 0.5):",
            font=("Segoe UI", 10, "bold"),
        ).pack(anchor="w", pady=(0, 4))

        rcols = ("step", "old_cap", "new_cap", "elems", "new_lf")
        self.tree_rehashes = ttk.Treeview(right_frame, columns=rcols, show="headings", height=5)
        self.tree_rehashes.heading("step", text="Rehash")
        self.tree_rehashes.heading("old_cap", text="Tam. Anterior")
        self.tree_rehashes.heading("new_cap", text="Nuevo Primo m")
        self.tree_rehashes.heading("elems", text="Elementos")
        self.tree_rehashes.heading("new_lf", text="Nuevo Load Factor")
        for col_id, w in zip(rcols, (65, 95, 105, 85, 120)):
            self.tree_rehashes.column(col_id, width=w, anchor="center")
        self.tree_rehashes.pack(fill="x", pady=(0, 10))

        ttk.Label(
            right_frame,
            text="Secuencia Completa Recuperada de Tabla Hash 1:",
            font=("Segoe UI", 10, "bold"),
        ).pack(anchor="w", pady=(0, 4))

        self.txt_sequence_view = tk.Text(
            right_frame,
            height=12,
            bg="#1e293b",
            fg="#e2e8f0",
            font=("Consolas", 10),
            wrap="char",
            relief="flat",
            padx=8,
            pady=8,
        )
        self.txt_sequence_view.pack(fill="both", expand=True)

        btn_row = ttk.Frame(right_frame)
        btn_row.pack(fill="x", pady=(8, 0))
        ttk.Button(
            btn_row,
            text="Analizar k-mers de esta Secuencia →",
            style="Accent.TButton",
            command=self._jump_to_kmer_analysis,
        ).pack(side="left", padx=(0, 6))
        ttk.Button(
            btn_row,
            text="Construir Grafo de esta Secuencia →",
            style="Accent.TButton",
            command=self._jump_to_graph_analysis,
        ).pack(side="left")

    # =========================================================================
    # PESTAÑA 2: TABLA HASH 2 (DOUBLE HASHING 4^n) & MAX-HEAP
    # =========================================================================
    def _build_tab_ht2_heap(self):
        controls = ttk.Frame(self.tab_ht2_heap, style="Card.TFrame", padding=10)
        controls.pack(fill="x", pady=(0, 10))

        ttk.Label(controls, text="Secuencia (ID):", style="Card.TLabel").pack(side="left", padx=(4, 6))
        self.combo_kmer_seq = ttk.Combobox(controls, width=28, font=("Consolas", 10))
        self.combo_kmer_seq.pack(side="left", padx=(0, 14))

        ttk.Label(controls, text="Tamaño k-mer (n):", style="Card.TLabel").pack(side="left", padx=(0, 6))
        self.spin_k = ttk.Spinbox(controls, from_=1, to=8, width=5, font=("Consolas", 10))
        self.spin_k.set("3")
        self.spin_k.pack(side="left", padx=(0, 14))

        ttk.Label(controls, text="Top-N (Max-Heap):", style="Card.TLabel").pack(side="left", padx=(0, 6))
        self.spin_topn = ttk.Spinbox(controls, from_=1, to=100, width=6, font=("Consolas", 10))
        self.spin_topn.set("15")
        self.spin_topn.pack(side="left", padx=(0, 14))

        ttk.Button(
            controls,
            text="Ejecutar Double Hashing + Max-Heap",
            style="Accent.TButton",
            command=self._run_kmer_and_heap,
        ).pack(side="left", padx=(0, 8))

        metrics_row = ttk.Frame(self.tab_ht2_heap)
        metrics_row.pack(fill="x", pady=(0, 10))

        c1, self.lbl_ht2_cap = self._create_metric_card(metrics_row, "Bins Disponibles (4^n)")
        c2, self.lbl_ht2_unique = self._create_metric_card(metrics_row, "k-mers Únicos Insertados")
        c3, self.lbl_ht2_total = self._create_metric_card(metrics_row, "Total k-mers Leídos (Paso 1)")
        c4, self.lbl_ht2_lf = self._create_metric_card(metrics_row, "Load Factor (Tabla Hash 2)")
        c5, self.lbl_ht2_coll = self._create_metric_card(metrics_row, "Colisiones (Double Hashing)")

        for card in (c1, c2, c3, c4, c5):
            card.pack(side="left", fill="x", expand=True, padx=4)

        paned = ttk.PanedWindow(self.tab_ht2_heap, orient="horizontal")
        paned.pack(fill="both", expand=True)

        # Izquierda: Tabla Hash 2 + Buscador de k-mer puntual + Max-Heap Top-N
        left_pane = ttk.Frame(paned, padding=(0, 0, 8, 0))
        paned.add(left_pane, weight=3)

        ksearch_frame = ttk.Frame(left_pane)
        ksearch_frame.pack(fill="x", pady=(0, 8))
        ttk.Label(ksearch_frame, text="Buscar k-mer en Tabla Hash 2:").pack(side="left", padx=(0, 6))
        self.entry_search_kmer = ttk.Entry(ksearch_frame, font=("Consolas", 10), width=12)
        self.entry_search_kmer.pack(side="left", padx=(0, 6))
        self.entry_search_kmer.insert(0, "ATG")
        self.entry_search_kmer.bind("<Return>", lambda _e: self._on_search_single_kmer())
        ttk.Button(
            ksearch_frame,
            text="Consultar Frecuencia",
            style="Accent.TButton",
            command=self._on_search_single_kmer,
        ).pack(side="left", padx=(0, 10))
        self.lbl_kmer_search_res = ttk.Label(ksearch_frame, text="", foreground="#38bdf8")
        self.lbl_kmer_search_res.pack(side="left")

        ttk.Label(
            left_pane,
            text="c) Ranking de los N Patrones (k-mers) Más Frecuentes extraídos del Max-Heap:",
            font=("Segoe UI", 10, "bold"),
        ).pack(anchor="w", pady=(0, 4))

        hcols = ("rank", "kmer", "count", "pct")
        self.tree_heap = ttk.Treeview(left_pane, columns=hcols, show="headings", height=14)
        self.tree_heap.heading("rank", text="Rank (Max-Heap)")
        self.tree_heap.heading("kmer", text="Key (k-mer)")
        self.tree_heap.heading("count", text="Valor (Frecuencia)")
        self.tree_heap.heading("pct", text="% del Total")
        self.tree_heap.column("rank", width=110, anchor="center")
        self.tree_heap.column("kmer", width=120, anchor="center")
        self.tree_heap.column("count", width=130, anchor="center")
        self.tree_heap.column("pct", width=110, anchor="center")
        self.tree_heap.pack(fill="both", expand=True)

        # Derecha: Gráfico de Barras del Max-Heap + Casillas de la Tabla Hash 2
        right_pane = ttk.Frame(paned, padding=(8, 0, 0, 0))
        paned.add(right_pane, weight=4)

        ttk.Label(
            right_pane,
            text="Visualización de Frecuencias (Max-Heap Top-N):",
            font=("Segoe UI", 10, "bold"),
        ).pack(anchor="w", pady=(0, 4))

        self.canvas_heap = tk.Canvas(right_pane, bg="#1e293b", height=240, highlightthickness=0)
        self.canvas_heap.pack(fill="x", pady=(0, 10))

        ttk.Label(
            right_pane,
            text="Inspección de Bins en Tabla Hash 2 (Double Hashing):",
            font=("Segoe UI", 10, "bold"),
        ).pack(anchor="w", pady=(0, 4))

        scols = ("slot", "state", "kmer", "count", "probes")
        self.tree_ht2_slots = ttk.Treeview(right_pane, columns=scols, show="headings", height=8)
        self.tree_ht2_slots.heading("slot", text="Bin #")
        self.tree_ht2_slots.heading("state", text="Estado")
        self.tree_ht2_slots.heading("kmer", text="Key (k-mer)")
        self.tree_ht2_slots.heading("count", text="Conteo")
        self.tree_ht2_slots.heading("probes", text="Sondeos (Double Hash)")
        for cid, w in zip(scols, (70, 90, 100, 80, 140)):
            self.tree_ht2_slots.column(cid, width=w, anchor="center")
        self.tree_ht2_slots.pack(fill="both", expand=True)

    # =========================================================================
    # PESTAÑA 3: GRAFO DE SECUENCIAS DE ADN (KHAN ET AL., 2024)
    # =========================================================================
    def _build_tab_graph(self):
        top_bar = ttk.Frame(self.tab_graph, style="Card.TFrame", padding=10)
        top_bar.pack(fill="x", pady=(0, 10))

        ttk.Label(top_bar, text="Secuencia FASTA:", style="Card.TLabel").pack(side="left", padx=(4, 6))
        self.combo_graph_seq = ttk.Combobox(top_bar, width=25, font=("Consolas", 10))
        self.combo_graph_seq.pack(side="left", padx=(0, 10))

        ttk.Label(top_bar, text="Parámetro β:", style="Card.TLabel").pack(side="left", padx=(0, 4))
        self.entry_beta = ttk.Entry(top_bar, width=6, font=("Consolas", 10))
        self.entry_beta.insert(0, "1.0")
        self.entry_beta.pack(side="left", padx=(0, 10))

        ttk.Button(
            top_bar,
            text="Construir Grafo de Secuencia FASTA",
            style="Accent.TButton",
            command=self._run_graph_from_fasta,
        ).pack(side="left", padx=(0, 12))

        ttk.Button(
            top_bar,
            text="★ Probar Ejemplo Exacto del PDF (TCATTAGCAGTTAGCAGC, β=1)",
            style="Accent.TButton",
            command=self._run_pdf_example_graph,
        ).pack(side="left")

        custom_bar = ttk.Frame(self.tab_graph)
        custom_bar.pack(fill="x", pady=(0, 8))
        ttk.Label(custom_bar, text="O ingresar cadena de ADN personalizada:").pack(side="left", padx=(0, 6))
        self.entry_custom_dna = ttk.Entry(custom_bar, font=("Consolas", 10), width=42)
        self.entry_custom_dna.pack(side="left", padx=(0, 8))
        self.entry_custom_dna.insert(0, "TCATTAGCAGTTAGCAGC")
        ttk.Button(
            custom_bar,
            text="Construir desde Cadena Personalizada",
            style="Accent.TButton",
            command=self._run_graph_from_custom,
        ).pack(side="left", padx=(0, 16))

        self.var_matrix_mode = tk.StringVar(value="active")
        ttk.Radiobutton(
            custom_bar,
            text="Matriz Compacta (Vértices presentes)",
            variable=self.var_matrix_mode,
            value="active",
            command=self._refresh_matrix_display,
        ).pack(side="left", padx=(0, 8))
        ttk.Radiobutton(
            custom_bar,
            text="Matriz Canónica 16×16 (AA..TT)",
            variable=self.var_matrix_mode,
            value="canonical",
            command=self._refresh_matrix_display,
        ).pack(side="left")

        self.lbl_graph_summary = ttk.Label(
            self.tab_graph,
            text="",
            foreground="#38bdf8",
            font=("Segoe UI", 10, "bold"),
        )
        self.lbl_graph_summary.pack(anchor="w", pady=(0, 6))

        paned = ttk.PanedWindow(self.tab_graph, orient="horizontal")
        paned.pack(fill="both", expand=True)

        # Izquierda: Tabla de posiciones de cada par (2-mer)
        left_frame = ttk.Frame(paned, padding=(0, 0, 8, 0))
        paned.add(left_frame, weight=1)

        ttk.Label(
            left_frame,
            text="Tabla de Pares (Posición → 2-mer):",
            font=("Segoe UI", 10, "bold"),
        ).pack(anchor="w", pady=(0, 4))

        pcols = ("pos", "pair")
        self.tree_pairs = ttk.Treeview(left_frame, columns=pcols, show="headings", height=16)
        self.tree_pairs.heading("pos", text="Posición de cada par")
        self.tree_pairs.heading("pair", text="2-mer")
        self.tree_pairs.column("pos", width=130, anchor="center")
        self.tree_pairs.column("pair", width=90, anchor="center")
        self.tree_pairs.pack(fill="both", expand=True)

        # Centro: Matriz de Adyacencia Ponderada
        mid_frame = ttk.Frame(paned, padding=(4, 0, 4, 0))
        paned.add(mid_frame, weight=3)

        ttk.Label(
            mid_frame,
            text="Matriz de Adyacencia Ponderada A[u, v]  (d = 1 / (j - i)^β):",
            font=("Segoe UI", 10, "bold"),
        ).pack(anchor="w", pady=(0, 4))

        self.txt_matrix = tk.Text(
            mid_frame,
            bg="#1e293b",
            fg="#f8fafc",
            font=("Consolas", 10),
            wrap="none",
            relief="flat",
            padx=10,
            pady=10,
        )
        scroll_x = ttk.Scrollbar(mid_frame, orient="horizontal", command=self.txt_matrix.xview)
        scroll_y = ttk.Scrollbar(mid_frame, orient="vertical", command=self.txt_matrix.yview)
        self.txt_matrix.configure(xscrollcommand=scroll_x.set, yscrollcommand=scroll_y.set)
        self.txt_matrix.pack(fill="both", expand=True)
        scroll_x.pack(fill="x")

        # Derecha: Dibujo del Grafo Dirigido en Canvas
        right_frame = ttk.Frame(paned, padding=(8, 0, 0, 0))
        paned.add(right_frame, weight=2)

        ttk.Label(
            right_frame,
            text="Grafo Dirigido de 2-mers (Loop Digraph G):",
            font=("Segoe UI", 10, "bold"),
        ).pack(anchor="w", pady=(0, 4))

        self.canvas_graph = tk.Canvas(right_frame, bg="#1e293b", width=340, highlightthickness=0)
        self.canvas_graph.pack(fill="both", expand=True)
        self.canvas_graph.bind("<Configure>", lambda _e: self._draw_graph_canvas())

    # =========================================================================
    # LOGICA DE INTERACCION CON EL MOTOR C++
    # =========================================================================
    def _refresh_initial_data(self):
        status = self.engine.send_command("STATUS")
        if not status.get("ok"):
            messagebox.showwarning("Aviso", status.get("message", "Error consultando estado"))
            return

        self.current_status = status
        self.sequences_meta = status.get("sequences", [])

        self.lbl_ht1_count.config(text=str(status.get("fastaCount", 0)))
        self.lbl_ht1_cap.config(text=str(status.get("ht1Capacity", 10)))
        self.lbl_ht1_lf.config(text=f"{status.get('ht1LoadFactor', 0.0):.4f}")
        self.lbl_ht1_rehashes.config(text=str(len(status.get("ht1Rehashes", []))))
        self.lbl_ht1_collisions.config(text=str(status.get("ht1Collisions", 0)))

        # Poblar tabla de secuencias
        for item in self.tree_sequences.get_children():
            self.tree_sequences.delete(item)
        seq_ids = []
        for s in self.sequences_meta:
            self.tree_sequences.insert(
                "",
                "end",
                values=(s["index"] + 1, s["id"], s["length"], s["preview"]),
            )
            seq_ids.append(s["id"])

        # Poblar historial de rehashes
        for item in self.tree_rehashes.get_children():
            self.tree_rehashes.delete(item)
        for idx, rh in enumerate(status.get("ht1Rehashes", []), start=1):
            self.tree_rehashes.insert(
                "",
                "end",
                values=(
                    f"#{idx}",
                    rh["oldCap"],
                    rh["newCap"],
                    rh["elements"],
                    f"{rh['loadFactor']:.4f}",
                ),
            )

        if seq_ids:
            self.combo_kmer_seq["values"] = seq_ids
            self.combo_kmer_seq.set(seq_ids[0])
            self.combo_graph_seq["values"] = seq_ids
            self.combo_graph_seq.set(seq_ids[0])
            self.entry_search_id.delete(0, "end")
            self.entry_search_id.insert(0, seq_ids[0])
            self._on_search_ht1()
            self._run_kmer_and_heap()

        self._run_pdf_example_graph()

    def _on_load_fasta_dialog(self):
        filepath = filedialog.askopenfilename(
            parent=self,
            title="Seleccionar archivo FASTA",
            initialdir=self.project_root,
            filetypes=[("Archivos FASTA", "*.fasta *.fa *.txt"), ("Todos los archivos", "*.*")],
        )
        if not filepath:
            return
        resp = self.engine.send_command(f"LOAD_FASTA {filepath}")
        if not resp.get("ok"):
            messagebox.showerror("Error", resp.get("message", "No se pudo cargar el archivo FASTA"))
            return
        self._refresh_initial_data()

    def _on_search_ht1(self):
        seq_id = self.entry_search_id.get().strip()
        if not seq_id:
            return
        resp = self.engine.send_command(f"SEARCH_SEQ {seq_id}")
        if resp.get("found"):
            seq_str = resp.get("sequence", "")
            probes = resp.get("comparisons", 1)
            a_cnt = seq_str.count("A")
            t_cnt = seq_str.count("T")
            c_cnt = seq_str.count("C")
            g_cnt = seq_str.count("G")
            self.lbl_search_result.config(
                text=f"✔ Encontrado en Tabla Hash 1: '{seq_id}' | Longitud: {len(seq_str)} pb | "
                f"Sondeos (Linear Probing): {probes} | A:{a_cnt} T:{t_cnt} C:{c_cnt} G:{g_cnt}",
                foreground="#4ade80",
            )
            self.txt_sequence_view.delete("1.0", "end")
            self.txt_sequence_view.insert("1.0", f">{seq_id}\n{seq_str}")
        else:
            probes = resp.get("comparisons", 0)
            self.lbl_search_result.config(
                text=f"✖ ID '{seq_id}' no encontrado en Tabla Hash 1 tras {probes} sondeo(s).",
                foreground="#f87171",
            )

    def _on_select_sequence_tree(self, _event):
        sel = self.tree_sequences.selection()
        if not sel:
            return
        vals = self.tree_sequences.item(sel[0], "values")
        if vals and len(vals) >= 2:
            seq_id = vals[1]
            self.entry_search_id.delete(0, "end")
            self.entry_search_id.insert(0, seq_id)
            self._on_search_ht1()

    def _jump_to_kmer_analysis(self):
        seq_id = self.entry_search_id.get().strip()
        if seq_id:
            self.combo_kmer_seq.set(seq_id)
            self.notebook.select(self.tab_ht2_heap)
            self._run_kmer_and_heap()

    def _jump_to_graph_analysis(self):
        seq_id = self.entry_search_id.get().strip()
        if seq_id:
            self.combo_graph_seq.set(seq_id)
            self.notebook.select(self.tab_graph)
            self._run_graph_from_fasta()

    def _run_kmer_and_heap(self):
        seq_id = self.combo_kmer_seq.get().strip()
        try:
            k = int(self.spin_k.get())
        except ValueError:
            k = 3
        try:
            top_n = int(self.spin_topn.get())
        except ValueError:
            top_n = 15

        resp = self.engine.send_command(f"ANALYZE_KMERS {seq_id} {k} {top_n}")
        if not resp.get("ok"):
            messagebox.showerror("Error", resp.get("message", "No se pudo analizar k-mers"))
            return

        cap = resp.get("ht2Capacity", 0)
        uniq = resp.get("ht2Unique", 0)
        tot = resp.get("ht2TotalValid", 0)
        lf = resp.get("ht2LoadFactor", 0.0)
        coll = resp.get("ht2Collisions", 0)

        self.lbl_ht2_cap.config(text=f"{cap} (4^{k})")
        self.lbl_ht2_unique.config(text=str(uniq))
        self.lbl_ht2_total.config(text=str(tot))
        self.lbl_ht2_lf.config(text=f"{lf:.4f}")
        self.lbl_ht2_coll.config(text=str(coll))

        for item in self.tree_heap.get_children():
            self.tree_heap.delete(item)

        top_kmers = resp.get("topKmers", [])
        for row in top_kmers:
            pct = (row["count"] / tot * 100.0) if tot > 0 else 0.0
            self.tree_heap.insert(
                "",
                "end",
                values=(f"#{row['rank']}", row["kmer"], row["count"], f"{pct:.2f}%"),
            )

        for item in self.tree_ht2_slots.get_children():
            self.tree_ht2_slots.delete(item)
        for slot in resp.get("slots", []):
            state_txt = "OCUPADO" if slot["occupied"] else "VACIO"
            self.tree_ht2_slots.insert(
                "",
                "end",
                values=(
                    slot["slot"],
                    state_txt,
                    slot["kmer"] if slot["occupied"] else "—",
                    slot["count"] if slot["occupied"] else 0,
                    slot["probes"] if slot["occupied"] else 0,
                ),
            )

        self._draw_heap_chart(top_kmers)
        if top_kmers:
            self.entry_search_kmer.delete(0, "end")
            self.entry_search_kmer.insert(0, top_kmers[0]["kmer"])
            self._on_search_single_kmer()

    def _on_search_single_kmer(self):
        kmer = self.entry_search_kmer.get().strip().upper()
        if not kmer:
            return
        resp = self.engine.send_command(f"SEARCH_KMER {kmer}")
        if not resp.get("ok"):
            return
        if resp.get("found"):
            self.lbl_kmer_search_res.config(
                text=f"✔ '{kmer}' ocurre {resp['count']} veces ({resp['comparisons']} sondeo(s) Double Hashing)",
                foreground="#4ade80",
            )
        else:
            self.lbl_kmer_search_res.config(
                text=f"✖ '{kmer}' no aparece (0 ocurrencias, {resp.get('comparisons', 0)} sondeos; k actual={resp.get('kLength')})",
                foreground="#f87171",
            )

    def _draw_heap_chart(self, top_kmers: list):
        self.canvas_heap.delete("all")
        if not top_kmers:
            return

        self.canvas_heap.update_idletasks()
        w = max(self.canvas_heap.winfo_width(), 480)
        h = max(self.canvas_heap.winfo_height(), 220)

        shown = top_kmers[:12]
        max_count = max(item["count"] for item in shown) if shown else 1
        pad_left, pad_right, pad_top, pad_bottom = 45, 20, 25, 42
        plot_w = w - pad_left - pad_right
        plot_h = h - pad_top - pad_bottom
        bar_gap = 8
        bar_w = max(18, (plot_w - bar_gap * len(shown)) / max(len(shown), 1))

        for i, item in enumerate(shown):
            x0 = pad_left + i * (bar_w + bar_gap)
            x1 = x0 + bar_w
            ratio = item["count"] / max_count if max_count > 0 else 0
            bar_h = ratio * plot_h
            y1 = h - pad_bottom
            y0 = y1 - bar_h

            color = "#38bdf8" if i == 0 else "#0284c7"
            self.canvas_heap.create_rectangle(x0, y0, x1, y1, fill=color, outline="")
            self.canvas_heap.create_text(
                (x0 + x1) / 2,
                y0 - 10,
                text=str(item["count"]),
                fill="#f8fafc",
                font=("Consolas", 9, "bold"),
            )
            self.canvas_heap.create_text(
                (x0 + x1) / 2,
                y1 + 16,
                text=item["kmer"],
                fill="#cbd5e1",
                font=("Consolas", 9),
            )

    def _get_beta_value(self) -> float:
        try:
            val = float(self.entry_beta.get().strip())
            return val if val > 0 else 1.0
        except ValueError:
            return 1.0

    def _run_graph_from_fasta(self):
        seq_id = self.combo_graph_seq.get().strip()
        if not seq_id:
            return
        beta = self._get_beta_value()
        resp = self.engine.send_command(f"BUILD_GRAPH_ID {seq_id} {beta}")
        if resp.get("ok"):
            self.current_graph_data = resp
            self._populate_graph_views()

    def _run_graph_from_custom(self):
        raw = self.entry_custom_dna.get().strip().upper()
        if not raw:
            return
        beta = self._get_beta_value()
        resp = self.engine.send_command(f"BUILD_GRAPH_RAW {raw} {beta}")
        if resp.get("ok"):
            self.current_graph_data = resp
            self._populate_graph_views()

    def _run_pdf_example_graph(self):
        self.entry_custom_dna.delete(0, "end")
        self.entry_custom_dna.insert(0, "TCATTAGCAGTTAGCAGC")
        self.entry_beta.delete(0, "end")
        self.entry_beta.insert(0, "1.0")
        self.var_matrix_mode.set("active")
        self._run_graph_from_custom()

    def _populate_graph_views(self):
        data = self.current_graph_data
        if not data:
            return

        for item in self.tree_pairs.get_children():
            self.tree_pairs.delete(item)
        for p in data.get("pairs", []):
            self.tree_pairs.insert("", "end", values=(p["pos"], p["pair"]))

        w_tc_at = data.get("w_TC_AT", 0.0)
        w_at_ca = data.get("w_AT_CA", 0.0)
        self.lbl_graph_summary.config(
            text=f"Origen: {data['label']}  |  Longitud: {data['seqLength']} pb  |  "
            f"Pares (2-mers): {data['pairCount']}  |  β = {data['beta']:.2f}  |  "
            f"w(TC, AT) = {w_tc_at:.5f}  |  w(AT, CA) = {w_at_ca:.5f}"
        )

        self._refresh_matrix_display()
        self._draw_graph_canvas()

    def _refresh_matrix_display(self):
        data = self.current_graph_data
        if not data:
            return

        mode = self.var_matrix_mode.get()
        if mode == "active":
            labels = data.get("activeLabels", [])
            matrix = data.get("activeMatrix", [])
            title = f"Matriz de Adyacencia de Vértices Activos ({len(labels)}×{len(labels)}, β={data['beta']}):\n\n"
        else:
            labels = data.get("canonicalLabels", [])
            matrix = data.get("canonicalMatrix", [])
            title = f"Matriz de Adyacencia Canónica Ω (16×16, β={data['beta']}):\n\n"

        lines = [title]
        header = "      " + "".join(f"{lbl:>9}" for lbl in labels)
        lines.append(header)
        lines.append("      " + "-" * (9 * len(labels)))

        for i, row_lbl in enumerate(labels):
            row_vals = "".join(
                f"{matrix[i][j]:>9.4f}" if matrix[i][j] > 0 else f"{'0':>9}"
                for j in range(len(labels))
            )
            lines.append(f"{row_lbl:>4} |{row_vals}")

        self.txt_matrix.delete("1.0", "end")
        self.txt_matrix.insert("1.0", "\n".join(lines))

    def _draw_graph_canvas(self):
        self.canvas_graph.delete("all")
        data = self.current_graph_data
        if not data:
            return

        labels = data.get("activeLabels", [])
        matrix = data.get("activeMatrix", [])
        pairs_list = data.get("pairs", [])
        beta = data.get("beta", 1.0)
        n = len(labels)
        if n == 0:
            return

        self.canvas_graph.update_idletasks()
        w = max(self.canvas_graph.winfo_width(), 340)
        h = max(self.canvas_graph.winfo_height(), 340)
        cx, cy = w / 2, h / 2
        radius = min(w, h) * 0.35
        node_r = 17

        label_to_idx = {lbl: i for i, lbl in enumerate(labels)}
        positions = {}
        # Si es el ejemplo del PDF (7 vertices con TA), colocar TA al centro como en la Figura 3 del paper
        if n == 7 and "TA" in label_to_idx:
            outer_labels = [lbl for lbl in labels if lbl != "TA"]
            positions[label_to_idx["TA"]] = (cx, cy)
            for k_idx, lbl in enumerate(outer_labels):
                angle = -math.pi / 2 + (2 * math.pi * k_idx / len(outer_labels))
                positions[label_to_idx[lbl]] = (
                    cx + radius * math.cos(angle),
                    cy + radius * math.sin(angle),
                )
        else:
            for i, lbl in enumerate(labels):
                angle = -math.pi / 2 + (2 * math.pi * i / n)
                x = cx + radius * math.cos(angle)
                y = cy + radius * math.sin(angle)
                positions[i] = (x, y)

        # Si hay pocos pares (<= 20, ej. TCATTAGCAGTTAGCAGC), dibujar el Loop Digraph de la Figura 3
        # con aristas paralelas curvas individuales para cada par (i, j)
        if 0 < len(pairs_list) <= 20:
            pair_occurrences = {}
            for idx_i in range(len(pairs_list)):
                for idx_j in range(idx_i + 1, len(pairs_list)):
                    u_lbl = pairs_list[idx_i]["pair"]
                    v_lbl = pairs_list[idx_j]["pair"]
                    pos_i = pairs_list[idx_i]["pos"]
                    pos_j = pairs_list[idx_j]["pos"]
                    if u_lbl not in label_to_idx or v_lbl not in label_to_idx:
                        continue
                    u = label_to_idx[u_lbl]
                    v = label_to_idx[v_lbl]
                    w_arc = 1.0 / ((pos_j - pos_i) ** beta)
                    key = (u, v)
                    pair_occurrences.setdefault(key, []).append(w_arc)

            for (u, v), weights in pair_occurrences.items():
                x1, y1 = positions[u]
                x2, y2 = positions[v]
                if u == v:
                    # Bucle (Self-loop) hacia afuera como en la Figura 3
                    dx, dy = x1 - cx, y1 - cy
                    d = math.hypot(dx, dy) or 1.0
                    ux, uy = dx / d, dy / d
                    lx = x1 + ux * (node_r + 14)
                    ly = y1 + uy * (node_r + 14)
                    self.canvas_graph.create_oval(
                        lx - 13, ly - 13, lx + 13, ly + 13, outline="#38bdf8", width=1.6
                    )
                    self.canvas_graph.create_text(
                        lx + ux * 16,
                        ly + uy * 16,
                        text=f"{sum(weights):.3g}",
                        fill="#38bdf8",
                        font=("Consolas", 8, "bold"),
                    )
                    continue

                dx, dy = x2 - x1, y2 - y1
                dist = math.hypot(dx, dy)
                if dist < 1e-3:
                    continue
                ux, uy = dx / dist, dy / dist
                nx, ny = -uy, ux  # Vector normal para curvar aristas paralelas
                sx, sy = x1 + ux * node_r, y1 + uy * node_r
                ex, ey = x2 - ux * node_r, y2 - uy * node_r

                m_count = len(weights)
                for k_arc, w_val in enumerate(weights):
                    offset = (k_arc - (m_count - 1) / 2.0) * 22.0
                    if (v, u) in pair_occurrences and offset == 0:
                        offset = 14.0
                    ctrl_x = (sx + ex) / 2 + nx * offset
                    ctrl_y = (sy + ey) / 2 + ny * offset
                    self.canvas_graph.create_line(
                        sx,
                        sy,
                        ctrl_x,
                        ctrl_y,
                        ex,
                        ey,
                        smooth=True,
                        fill="#64748b",
                        width=1.4,
                        arrow=tk.LAST,
                        arrowshape=(8, 10, 4),
                    )
                    self.canvas_graph.create_text(
                        ctrl_x,
                        ctrl_y - 6,
                        text=f"{w_val:.4g}",
                        fill="#cbd5e1",
                        font=("Consolas", 8),
                    )
        else:
            # Para secuencias largas (HIV-2), mostrar las 35 aristas mas pesadas
            edges = []
            for i in range(n):
                for j in range(n):
                    if matrix[i][j] > 0:
                        edges.append((matrix[i][j], i, j))
            edges.sort(reverse=True)

            for weight, u, v in edges[:35]:
                x1, y1 = positions[u]
                x2, y2 = positions[v]
                if u == v:
                    self.canvas_graph.create_oval(
                        x1 - 14, y1 - 32, x1 + 14, y1 - 6, outline="#38bdf8", width=1.5
                    )
                    continue
                dx, dy = x2 - x1, y2 - y1
                dist = math.hypot(dx, dy)
                if dist < 1e-3:
                    continue
                ux, uy = dx / dist, dy / dist
                sx, sy = x1 + ux * node_r, y1 + uy * node_r
                ex, ey = x2 - ux * node_r, y2 - uy * node_r
                self.canvas_graph.create_line(
                    sx,
                    sy,
                    ex,
                    ey,
                    fill="#475569",
                    width=1.4,
                    arrow=tk.LAST,
                    arrowshape=(8, 10, 4),
                )

        # Dibujar vertices
        for i, lbl in enumerate(labels):
            x, y = positions[i]
            self.canvas_graph.create_oval(
                x - node_r,
                y - node_r,
                x + node_r,
                y + node_r,
                fill="#0284c7",
                outline="#38bdf8",
                width=2,
            )
            self.canvas_graph.create_text(
                x,
                y,
                text=lbl,
                fill="#ffffff",
                font=("Consolas", 10, "bold"),
            )

    def _on_close(self):
        if hasattr(self, "engine") and self.engine:
            self.engine.close()
        self.destroy()


if __name__ == "__main__":
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    app = GenomicAnalyzerApp(root_dir)
    app.mainloop()
