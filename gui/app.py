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
            errors="replace",
            bufsize=1,
            creationflags=creationflags,
        )
        if self.proc.poll() is not None:
            err = ""
            try:
                err = self.proc.stderr.read().strip()
            except Exception:
                pass
            raise RuntimeError(
                f"El motor C++ terminó inesperadamente (código {self.proc.returncode})."
                + (f"\nDetalle: {err}" if err else "")
            )

    def send_command(self, command_line: str) -> dict:
        if self.proc is None or self.proc.poll() is not None:
            self._start_engine()
        try:
            self.proc.stdin.write(command_line.strip() + "\n")
            self.proc.stdin.flush()
            response_line = self.proc.stdout.readline()
        except Exception as exc:
            return {"ok": False, "message": f"Error de comunicación con el motor C++: {exc}"}
        if not response_line:
            rc = self.proc.poll()
            err = ""
            if rc is not None and self.proc.stderr:
                try:
                    err = self.proc.stderr.read().strip()
                except Exception:
                    pass
            msg = "Sin respuesta del motor C++"
            if rc is not None:
                msg += f" (salió con código {rc})"
            if err:
                msg += f": {err}"
            return {"ok": False, "message": msg}
        try:
            return json.loads(response_line)
        except json.JSONDecodeError as exc:
            return {"ok": False, "message": f"Error decodificando JSON: {exc}"}

    def close(self):
        if self.proc and self.proc.poll() is None:
            try:
                if self.proc.stdin:
                    self.proc.stdin.write("QUIT\n")
                    self.proc.stdin.flush()
                self.proc.wait(timeout=1.0)
            except subprocess.TimeoutExpired:
                self.proc.terminate()
            except Exception:
                try:
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
        self.current_top_kmers = []
        self.current_k = 3

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
        bg_input = "#0f172a"
        border_col = "#475569"
        fg_main = "#f8fafc"
        accent = "#38bdf8"
        accent_active = "#0284c7"

        self.option_add("*TCombobox*Listbox.background", bg_card)
        self.option_add("*TCombobox*Listbox.foreground", fg_main)
        self.option_add("*TCombobox*Listbox.selectBackground", accent_active)
        self.option_add("*TCombobox*Listbox.selectForeground", "#ffffff")

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
            "TEntry",
            fieldbackground=bg_input,
            foreground=fg_main,
            insertcolor=fg_main,
            bordercolor=border_col,
            lightcolor=border_col,
            darkcolor=border_col,
            padding=(6, 4),
        )
        style.map(
            "TEntry",
            fieldbackground=[("focus", bg_input), ("readonly", bg_card)],
            bordercolor=[("focus", accent)],
        )

        style.configure(
            "TCombobox",
            fieldbackground=bg_input,
            background="#334155",
            foreground=fg_main,
            arrowcolor=fg_main,
            bordercolor=border_col,
            lightcolor=border_col,
            darkcolor=border_col,
            selectbackground=accent_active,
            selectforeground="#ffffff",
            padding=(6, 4),
        )
        style.map(
            "TCombobox",
            fieldbackground=[("readonly", bg_input), ("focus", bg_input)],
            foreground=[("readonly", fg_main)],
            background=[("active", border_col)],
            bordercolor=[("focus", accent)],
        )

        style.configure(
            "TSpinbox",
            fieldbackground=bg_input,
            background="#334155",
            foreground=fg_main,
            arrowcolor=fg_main,
            insertcolor=fg_main,
            bordercolor=border_col,
            lightcolor=border_col,
            darkcolor=border_col,
            padding=(4, 4),
        )
        style.map(
            "TSpinbox",
            fieldbackground=[("focus", bg_input)],
            background=[("active", border_col)],
            bordercolor=[("focus", accent)],
        )

        style.configure(
            "TRadiobutton",
            background=bg_dark,
            foreground=fg_main,
            indicatorcolor=bg_card,
            indicatorrelief="flat",
        )
        style.map(
            "TRadiobutton",
            background=[("active", bg_dark)],
            foreground=[("active", accent)],
            indicatorcolor=[("selected", accent), ("active", "#334155")],
        )

        style.configure(
            "Vertical.TScrollbar",
            background="#334155",
            troughcolor=bg_dark,
            bordercolor=bg_card,
            arrowcolor=fg_main,
        )
        style.configure(
            "Horizontal.TScrollbar",
            background="#334155",
            troughcolor=bg_dark,
            bordercolor=bg_card,
            arrowcolor=fg_main,
        )
        style.map(
            "TScrollbar",
            background=[("active", border_col), ("pressed", accent_active)],
        )

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

    def _create_tree_container(self, parent, pack_fill="both", pack_expand=True, pack_pady=(0, 0)):
        container = ttk.Frame(parent)
        container.pack(fill=pack_fill, expand=pack_expand, pady=pack_pady)
        container.columnconfigure(0, weight=1)
        container.rowconfigure(0, weight=1)
        return container

    def _attach_tree_scrollbar(self, container, tree):
        vscroll = ttk.Scrollbar(container, orient="vertical", command=tree.yview)
        tree.configure(yscrollcommand=vscroll.set)
        tree.grid(row=0, column=0, sticky="nsew")
        vscroll.grid(row=0, column=1, sticky="ns")
        return vscroll

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
        seq_container = self._create_tree_container(left_frame)
        self.tree_sequences = ttk.Treeview(seq_container, columns=cols, show="headings", height=12)
        self.tree_sequences.heading("idx", text="#")
        self.tree_sequences.heading("id", text="ID Secuencia (Key)")
        self.tree_sequences.heading("length", text="Longitud (pb)")
        self.tree_sequences.heading("preview", text="Prefijo Cadena ADN (Valor)")
        self.tree_sequences.column("idx", width=45, anchor="center")
        self.tree_sequences.column("id", width=210, anchor="w")
        self.tree_sequences.column("length", width=95, anchor="center")
        self.tree_sequences.column("preview", width=260, anchor="w")
        self._attach_tree_scrollbar(seq_container, self.tree_sequences)
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
        rh_container = self._create_tree_container(right_frame, pack_fill="x", pack_expand=False, pack_pady=(0, 10))
        self.tree_rehashes = ttk.Treeview(rh_container, columns=rcols, show="headings", height=5)
        self.tree_rehashes.heading("step", text="Rehash")
        self.tree_rehashes.heading("old_cap", text="Tam. Anterior")
        self.tree_rehashes.heading("new_cap", text="Nuevo Primo m")
        self.tree_rehashes.heading("elems", text="Elementos")
        self.tree_rehashes.heading("new_lf", text="Nuevo Load Factor")
        for col_id, w in zip(rcols, (65, 95, 105, 85, 120)):
            self.tree_rehashes.column(col_id, width=w, anchor="center")
        self._attach_tree_scrollbar(rh_container, self.tree_rehashes)

        ttk.Label(
            right_frame,
            text="Secuencia Completa Recuperada de Tabla Hash 1:",
            font=("Segoe UI", 10, "bold"),
        ).pack(anchor="w", pady=(0, 4))

        seq_txt_container = self._create_tree_container(right_frame)
        self.txt_sequence_view = tk.Text(
            seq_txt_container,
            height=12,
            bg="#1e293b",
            fg="#e2e8f0",
            font=("Consolas", 10),
            wrap="char",
            relief="flat",
            padx=8,
            pady=8,
        )
        seq_txt_scroll = ttk.Scrollbar(seq_txt_container, orient="vertical", command=self.txt_sequence_view.yview)
        self.txt_sequence_view.configure(yscrollcommand=seq_txt_scroll.set)
        self.txt_sequence_view.grid(row=0, column=0, sticky="nsew")
        seq_txt_scroll.grid(row=0, column=1, sticky="ns")

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
        self.combo_kmer_seq = ttk.Combobox(controls, width=24, font=("Consolas", 10))
        self.combo_kmer_seq.pack(side="left", padx=(0, 10))

        ttk.Label(controls, text="Tamaño k-mer (n):", style="Card.TLabel").pack(side="left", padx=(0, 6))
        self.spin_k = ttk.Spinbox(controls, from_=1, to=10, width=4, font=("Consolas", 10))
        self.spin_k.set("3")
        self.spin_k.pack(side="left", padx=(0, 10))

        ttk.Label(controls, text="Top-N (Max-Heap):", style="Card.TLabel").pack(side="left", padx=(0, 6))
        self.spin_topn = ttk.Spinbox(controls, from_=1, to=100, width=5, font=("Consolas", 10))
        self.spin_topn.set("15")
        self.spin_topn.pack(side="left", padx=(0, 10))

        ttk.Button(
            controls,
            text="Ejecutar Double Hashing + Max-Heap",
            style="Accent.TButton",
            command=self._run_kmer_and_heap,
        ).pack(side="left", padx=(0, 8))

        ttk.Button(
            controls,
            text="★ Ejemplo PDF (ATGGGATGTCTTG)",
            style="Accent.TButton",
            command=self._run_pdf_example_kmers,
        ).pack(side="left")

        metrics_row = ttk.Frame(self.tab_ht2_heap)
        metrics_row.pack(fill="x", pady=(0, 6))

        c1, self.lbl_ht2_cap = self._create_metric_card(metrics_row, "Bins Disponibles (4^n)")
        c2, self.lbl_ht2_unique = self._create_metric_card(metrics_row, "k-mers Únicos Insertados")
        c3, self.lbl_ht2_total = self._create_metric_card(metrics_row, "Total k-mers Leídos (Paso 1)")
        c4, self.lbl_ht2_lf = self._create_metric_card(metrics_row, "Load Factor (Tabla Hash 2)")
        c5, self.lbl_ht2_coll = self._create_metric_card(metrics_row, "Colisiones (Double Hashing)")

        for card in (c1, c2, c3, c4, c5):
            card.pack(side="left", fill="x", expand=True, padx=4)

        self.lbl_ht2_contrast = ttk.Label(
            self.tab_ht2_heap,
            text="Contrastación DH vs Linear Probing y Max-Heap disponible tras ejecutar análisis.",
            foreground="#94a3b8",
            font=("Segoe UI", 9, "bold"),
        )
        self.lbl_ht2_contrast.pack(anchor="w", padx=4, pady=(0, 8))

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
        heap_container = self._create_tree_container(left_pane)
        self.tree_heap = ttk.Treeview(heap_container, columns=hcols, show="headings", height=14)
        self.tree_heap.heading("rank", text="Rank (Max-Heap)")
        self.tree_heap.heading("kmer", text="Key (k-mer)")
        self.tree_heap.heading("count", text="Valor (Frecuencia)")
        self.tree_heap.heading("pct", text="% del Total")
        self.tree_heap.column("rank", width=110, anchor="center")
        self.tree_heap.column("kmer", width=120, anchor="center")
        self.tree_heap.column("count", width=130, anchor="center")
        self.tree_heap.column("pct", width=110, anchor="center")
        self._attach_tree_scrollbar(heap_container, self.tree_heap)

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
        self.canvas_heap.bind("<Configure>", lambda _e: self._draw_heap_chart(self.current_top_kmers))

        ttk.Label(
            right_pane,
            text="Inspección de Bins en Tabla Hash 2 (Double Hashing):",
            font=("Segoe UI", 10, "bold"),
        ).pack(anchor="w", pady=(0, 4))

        scols = ("slot", "state", "kmer", "count", "probes")
        ht2_container = self._create_tree_container(right_pane)
        self.tree_ht2_slots = ttk.Treeview(ht2_container, columns=scols, show="headings", height=8)
        self.tree_ht2_slots.heading("slot", text="Bin #")
        self.tree_ht2_slots.heading("state", text="Estado")
        self.tree_ht2_slots.heading("kmer", text="Key (k-mer)")
        self.tree_ht2_slots.heading("count", text="Conteo")
        self.tree_ht2_slots.heading("probes", text="Sondeos (Double Hash)")
        for cid, w in zip(scols, (70, 90, 100, 80, 140)):
            self.tree_ht2_slots.column(cid, width=w, anchor="center")
        self._attach_tree_scrollbar(ht2_container, self.tree_ht2_slots)

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
        pairs_container = self._create_tree_container(left_frame)
        self.tree_pairs = ttk.Treeview(pairs_container, columns=pcols, show="headings", height=16)
        self.tree_pairs.heading("pos", text="Posición de cada par")
        self.tree_pairs.heading("pair", text="2-mer")
        self.tree_pairs.column("pos", width=130, anchor="center")
        self.tree_pairs.column("pair", width=90, anchor="center")
        self._attach_tree_scrollbar(pairs_container, self.tree_pairs)

        # Centro: Matriz de Adyacencia Ponderada
        mid_frame = ttk.Frame(paned, padding=(4, 0, 4, 0))
        paned.add(mid_frame, weight=3)

        ttk.Label(
            mid_frame,
            text="Matriz de Adyacencia Ponderada A[u, v]  (d = 1 / (j - i)^β):",
            font=("Segoe UI", 10, "bold"),
        ).pack(anchor="w", pady=(0, 4))

        matrix_container = self._create_tree_container(mid_frame)
        self.txt_matrix = tk.Text(
            matrix_container,
            bg="#1e293b",
            fg="#f8fafc",
            font=("Consolas", 10),
            wrap="none",
            relief="flat",
            padx=10,
            pady=10,
        )
        scroll_x = ttk.Scrollbar(matrix_container, orient="horizontal", command=self.txt_matrix.xview)
        scroll_y = ttk.Scrollbar(matrix_container, orient="vertical", command=self.txt_matrix.yview)
        self.txt_matrix.configure(xscrollcommand=scroll_x.set, yscrollcommand=scroll_y.set)
        self.txt_matrix.grid(row=0, column=0, sticky="nsew")
        scroll_y.grid(row=0, column=1, sticky="ns")
        scroll_x.grid(row=1, column=0, sticky="ew")

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
        self._apply_status_data(status)
        self._run_pdf_example_graph()

    def _apply_status_data(self, status: dict):
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
        else:
            self.combo_kmer_seq["values"] = []
            self.combo_graph_seq["values"] = []
            self.txt_sequence_view.delete("1.0", "end")
            self.lbl_search_result.config(
                text="No hay secuencias cargadas en Tabla Hash 1.",
                foreground="#94a3b8",
            )

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
        self._apply_status_data(resp)

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
            other_cnt = len(seq_str) - (a_cnt + t_cnt + c_cnt + g_cnt)
            base_summary = f"A:{a_cnt} T:{t_cnt} C:{c_cnt} G:{g_cnt}"
            if other_cnt > 0:
                base_summary += f" | IUPAC/otros:{other_cnt}"
            self.lbl_search_result.config(
                text=f"✔ Encontrado en Tabla Hash 1: '{seq_id}' | Longitud: {len(seq_str)} pb | "
                f"Sondeos (Linear Probing): {probes} | {base_summary}",
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
            self.txt_sequence_view.delete("1.0", "end")

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
        if not seq_id:
            messagebox.showwarning("Secuencia requerida", "Seleccione o ingrese el ID de una secuencia.")
            return
        try:
            k = int(self.spin_k.get().strip())
        except ValueError:
            messagebox.showwarning("Parámetro k inválido", "El tamaño k-mer debe ser un entero entre 1 y 10.")
            return
        if k < 1 or k > 10:
            messagebox.showwarning("Parámetro k inválido", "El tamaño k-mer debe estar entre 1 y 10.")
            return
        try:
            top_n = int(self.spin_topn.get().strip())
        except ValueError:
            messagebox.showwarning("Parámetro Top-N inválido", "La cantidad Top-N debe ser un entero >= 1.")
            return
        if top_n < 1:
            messagebox.showwarning("Parámetro Top-N inválido", "La cantidad Top-N debe ser mayor o igual a 1.")
            return

        resp = self.engine.send_command(f"ANALYZE_KMERS {seq_id} {k} {top_n}")
        if not resp.get("ok"):
            messagebox.showerror("Error", resp.get("message", "No se pudo analizar k-mers"))
            return
        self._apply_kmer_analysis_response(resp)

    def _run_pdf_example_kmers(self):
        self.spin_k.set("3")
        self.spin_topn.set("10")
        resp = self.engine.send_command("ANALYZE_KMERS_RAW ATGGGATGTCTTG 3 10")
        if not resp.get("ok"):
            messagebox.showerror("Error", resp.get("message", "No se pudo analizar el ejemplo del PDF"))
            return
        self._apply_kmer_analysis_response(resp)

    def _apply_kmer_analysis_response(self, resp: dict):
        k = int(resp.get("k", 3))
        cap = resp.get("ht2Capacity", 0)
        uniq = resp.get("ht2Unique", 0)
        tot = resp.get("ht2TotalValid", 0)
        lf = resp.get("ht2LoadFactor", 0.0)
        coll = resp.get("ht2Collisions", 0)
        dh_probes = resp.get("ht2TotalProbes", 0)
        dh_max = resp.get("ht2MaxProbes", 0)
        lp_probes = resp.get("lpTotalProbes", 0)
        lp_coll = resp.get("lpCollisions", 0)
        lp_max = resp.get("lpMaxProbes", 0)
        heap_cmp = resp.get("heapComparisons", 0)

        self.current_k = k
        self.lbl_ht2_cap.config(text=f"{cap} (4^{k})")
        self.lbl_ht2_unique.config(text=str(uniq))
        self.lbl_ht2_total.config(text=str(tot))
        self.lbl_ht2_lf.config(text=f"{lf:.4f}")
        self.lbl_ht2_coll.config(text=str(coll))

        dh_avg = (dh_probes / tot) if tot > 0 else 0.0
        lp_avg = (lp_probes / tot) if tot > 0 else 0.0
        self.lbl_ht2_contrast.config(
            text=(
                f"Contrastación ({resp.get('seqId', '')}, k={k}, ventanas={tot}): "
                f"Double Hashing = {dh_probes} sondas (prom {dh_avg:.4f}, col {coll}, max {dh_max})  |  "
                f"Linear Probing = {lp_probes} sondas (prom {lp_avg:.4f}, col {lp_coll}, max {lp_max})  |  "
                f"Max-Heap comparaciones = {heap_cmp}"
            ),
            foreground="#cbd5e1",
        )

        for item in self.tree_heap.get_children():
            self.tree_heap.delete(item)

        top_kmers = resp.get("topKmers", [])
        self.current_top_kmers = top_kmers
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
        else:
            self.lbl_kmer_search_res.config(text="")

    def _on_search_single_kmer(self):
        kmer = self.entry_search_kmer.get().strip().upper()
        if not kmer:
            messagebox.showwarning("k-mer requerido", "Ingrese un k-mer válido compuesto por A, T, C y G.")
            return
        if len(kmer) != self.current_k:
            messagebox.showwarning(
                "Longitud inválida",
                f"El k-mer '{kmer}' tiene longitud {len(kmer)}, pero la Tabla Hash 2 usa k={self.current_k}.",
            )
            return
        if any(c not in "ATCG" for c in kmer):
            messagebox.showwarning(
                "Caracteres inválidos",
                f"El k-mer '{kmer}' debe contener únicamente bases A, T, C y G.",
            )
            return

        resp = self.engine.send_command(f"SEARCH_KMER {kmer}")
        if not resp.get("ok"):
            messagebox.showerror("Error", resp.get("message", "No se pudo consultar el k-mer"))
            return
        if resp.get("found"):
            self.lbl_kmer_search_res.config(
                text=f"✔ '{kmer}' ocurre {resp['count']} veces ({resp['comparisons']} sondeo(s) Double Hashing)",
                foreground="#4ade80",
            )
        else:
            self.lbl_kmer_search_res.config(
                text=(
                    f"✖ '{kmer}' no aparece (0 ocurrencias — posible nulómero, "
                    f"{resp.get('comparisons', 0)} sondeos; k actual={resp.get('kLength')})"
                ),
                foreground="#f87171",
            )

    def _draw_heap_chart(self, top_kmers: list):
        self.canvas_heap.delete("all")
        if not top_kmers:
            return

        self.canvas_heap.update_idletasks()
        w = max(self.canvas_heap.winfo_width(), 480)
        h = max(self.canvas_heap.winfo_height(), 220)

        max_bars = 10 if self.current_k >= 6 else 12
        shown = top_kmers[:max_bars]
        max_count = max(item["count"] for item in shown) if shown else 1
        pad_left, pad_right, pad_top, pad_bottom = 45, 20, 25, (48 if self.current_k >= 6 else 40)
        plot_w = w - pad_left - pad_right
        plot_h = max(h - pad_top - pad_bottom, 40)
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
            if self.current_k >= 6:
                lbl_y = y1 + (13 if i % 2 == 0 else 29)
                lbl_font = ("Consolas", 8)
            else:
                lbl_y = y1 + 16
                lbl_font = ("Consolas", 9)
            self.canvas_heap.create_text(
                (x0 + x1) / 2,
                lbl_y,
                text=item["kmer"],
                fill="#cbd5e1",
                font=lbl_font,
            )

    def _get_beta_value(self):
        raw = self.entry_beta.get().strip()
        try:
            val = float(raw)
        except ValueError:
            messagebox.showwarning(
                "Parámetro β inválido",
                f"El parámetro β ('{raw}') debe ser un número real positivo (> 0).",
            )
            return None
        if not math.isfinite(val) or val <= 0.0:
            messagebox.showwarning(
                "Parámetro β inválido",
                f"El parámetro β ({val}) debe ser finito y estrictamente mayor que 0.",
            )
            return None
        return val

    def _run_graph_from_fasta(self):
        seq_id = self.combo_graph_seq.get().strip()
        if not seq_id:
            messagebox.showwarning("Secuencia requerida", "Seleccione una secuencia FASTA válida.")
            return
        beta = self._get_beta_value()
        if beta is None:
            return
        resp = self.engine.send_command(f"BUILD_GRAPH_ID {seq_id} {beta}")
        if not resp.get("ok"):
            messagebox.showerror("Error al construir grafo", resp.get("message", "No se pudo construir el grafo"))
            return
        self.current_graph_data = resp
        self._populate_graph_views()

    def _run_graph_from_custom(self):
        raw = "".join(self.entry_custom_dna.get().split()).upper()
        if len(raw) < 2:
            messagebox.showwarning(
                "Cadena insuficiente",
                "Ingrese al menos 2 nucleótidos válidos (A, T, C, G) para formar dinucleótidos.",
            )
            return
        beta = self._get_beta_value()
        if beta is None:
            return
        resp = self.engine.send_command(f"BUILD_GRAPH_RAW {raw} {beta}")
        if not resp.get("ok"):
            messagebox.showerror("Error al construir grafo", resp.get("message", "No se pudo construir el grafo"))
            return
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
        if labels and matrix:
            header = "      " + "".join(f" {lbl:>9}" for lbl in labels)
            lines.append(header)
            lines.append("      " + "-" * (10 * len(labels)))

            for i, row_lbl in enumerate(labels):
                row_vals = "".join(
                    f" {matrix[i][j]:>9.4f}" if matrix[i][j] > 0 else f" {'0':>9}"
                    for j in range(len(labels))
                )
                lines.append(f"{row_lbl:>4} |{row_vals}")
        else:
            lines.append("(Sin matriz disponible)")

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
        total_pairs = data.get("pairCount", len(pairs_list))
        beta = data.get("beta", 1.0)
        n = len(labels)
        if n == 0:
            return

        self.canvas_graph.update_idletasks()
        w = max(self.canvas_graph.winfo_width(), 340)
        h = max(self.canvas_graph.winfo_height(), 340)
        cx, cy = w / 2, h / 2 - 6
        radius = min(w, h) * 0.33
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
                    dx, dy = x1 - cx, y1 - cy
                    d = math.hypot(dx, dy) or 1.0
                    ux, uy = dx / d, dy / d
                    lx = x1 + ux * (node_r + 12)
                    ly = y1 + uy * (node_r + 12)
                    self.canvas_graph.create_oval(
                        lx - 12, ly - 12, lx + 12, ly + 12, outline="#38bdf8", width=1.5
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

        if 0 < total_pairs <= 20 and len(pairs_list) <= 20:
            mode_txt = f"Modo: Loop Digraph exacto ({total_pairs} pares)"
        else:
            mode_txt = f"Modo: Subgrafo Top-35 aristas ({total_pairs} pares totales)"
        self.canvas_graph.create_text(
            10,
            h - 12,
            anchor="w",
            text=mode_txt,
            fill="#94a3b8",
            font=("Segoe UI", 8),
        )

    def _on_close(self):
        if hasattr(self, "engine") and self.engine:
            self.engine.close()
        self.destroy()


if __name__ == "__main__":
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    app = GenomicAnalyzerApp(root_dir)
    app.mainloop()
