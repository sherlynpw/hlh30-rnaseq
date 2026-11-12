#!/usr/bin/env bash
# Extend gene_db.sqlite with WormCat and GO tables (+ FTS5).
# Idempotent — drops and recreates the tables it owns; leaves annotations / entrez tables untouched.
#
# Usage:
#   scripts/build_extended_tables.sh [DB] [WORMCAT_CSV] [CACHE_DIR]

set -euo pipefail

DB="${1:-gene_db.sqlite}"
WORMCAT_CSV="${2:-whole_genome_v2_nov-11-2021.csv}"
GENEIDS_TXT_GZ="${GENEIDS_TXT_GZ:-c_elegans.PRJNA13758.WS295.geneIDs.txt.gz}"
CACHE="${3:-.cache/go}"
ORTHO_CACHE="${ORTHO_CACHE:-.cache/orthology}"
BIOGRID_VERSION="${BIOGRID_VERSION:-5.0.257}"
BIOGRID_CACHE="${BIOGRID_CACHE:-.cache/biogrid}"

GAF_URL="https://current.geneontology.org/annotations/wb.gaf.gz"
OBO_URL="https://current.geneontology.org/ontology/go.obo"
BIOGRID_URL="https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-${BIOGRID_VERSION}/BIOGRID-ORGANISM-${BIOGRID_VERSION}.tab3.zip"
ALLIANCE_ORTHO_URL="https://fms.alliancegenome.org/download/ORTHOLOGY-ALLIANCE_COMBINED.tsv.gz"
ALLIANCE_PHENO_URL="https://fms.alliancegenome.org/download/PHENOTYPE_WB.json.gz"
WBPHENO_OBO_URL="https://raw.githubusercontent.com/obophenotype/c-elegans-phenotype-ontology/master/wbphenotype.obo"
PHENO_CACHE="${PHENO_CACHE:-.cache/phenotype}"

# WormBase anatomy ontology (WBbt) + gene-anatomy associations for tissue
# enrichment (TEA-style). downloads.wormbase.org is bot-blocked by Cloudflare,
# so we pull from publicly accessible mirrors:
#   * WBbt OBO from the obophenotype GitHub repo
#   * gene-anatomy expression rows from Alliance of Genome Resources
ANATOMY_OBO_URL="https://raw.githubusercontent.com/obophenotype/c-elegans-gross-anatomy-ontology/master/wbbt.obo"
ANATOMY_ASSOC_URL="https://fms.alliancegenome.org/download/EXPRESSION-ALLIANCE_WB.tsv"
ANATOMY_CACHE="${ANATOMY_CACHE:-.cache/anatomy}"

mkdir -p "$CACHE"

if [[ ! -f "$CACHE/wb.gaf.gz" ]]; then
  echo "Downloading $GAF_URL ..." >&2
  curl -fsSL "$GAF_URL" -o "$CACHE/wb.gaf.gz"
fi

if [[ ! -f "$CACHE/go.obo" ]]; then
  echo "Downloading $OBO_URL ..." >&2
  curl -fsSL "$OBO_URL" -o "$CACHE/go.obo"
fi

mkdir -p "$ORTHO_CACHE"
ALLIANCE_ORTHO_TSV_GZ="$ORTHO_CACHE/ORTHOLOGY-ALLIANCE_COMBINED.tsv.gz"
if [[ ! -f "$ALLIANCE_ORTHO_TSV_GZ" ]]; then
  echo "Downloading $ALLIANCE_ORTHO_URL ..." >&2
  curl -fsSL "$ALLIANCE_ORTHO_URL" -o "$ALLIANCE_ORTHO_TSV_GZ"
fi

mkdir -p "$PHENO_CACHE"
ALLIANCE_PHENO_JSON_GZ="$PHENO_CACHE/PHENOTYPE_WB.json.gz"
WBPHENO_OBO="$PHENO_CACHE/wbphenotype.obo"
if [[ ! -f "$ALLIANCE_PHENO_JSON_GZ" ]]; then
  echo "Downloading $ALLIANCE_PHENO_URL ..." >&2
  curl -fsSL "$ALLIANCE_PHENO_URL" -o "$ALLIANCE_PHENO_JSON_GZ"
fi
if [[ ! -f "$WBPHENO_OBO" ]]; then
  echo "Downloading $WBPHENO_OBO_URL ..." >&2
  curl -fsSL "$WBPHENO_OBO_URL" -o "$WBPHENO_OBO"
fi

mkdir -p "$ANATOMY_CACHE"
ANATOMY_OBO="$ANATOMY_CACHE/wbbt.obo"
ANATOMY_ASSOC="$ANATOMY_CACHE/EXPRESSION-ALLIANCE_WB.tsv"
if [[ ! -f "$ANATOMY_OBO" ]]; then
  echo "Downloading $ANATOMY_OBO_URL ..." >&2
  curl -fsSL "$ANATOMY_OBO_URL" -o "$ANATOMY_OBO"
fi
if [[ ! -f "$ANATOMY_ASSOC" ]]; then
  echo "Downloading $ANATOMY_ASSOC_URL ..." >&2
  curl -fsSL "$ANATOMY_ASSOC_URL" -o "$ANATOMY_ASSOC"
fi

mkdir -p "$BIOGRID_CACHE"
BIOGRID_TAB="$BIOGRID_CACHE/BIOGRID-ORGANISM-Caenorhabditis_elegans-${BIOGRID_VERSION}.tab3.txt"
if [[ ! -f "$BIOGRID_TAB" ]]; then
  BIOGRID_ZIP="$BIOGRID_CACHE/BIOGRID-ORGANISM-${BIOGRID_VERSION}.tab3.zip"
  if [[ ! -f "$BIOGRID_ZIP" ]]; then
    echo "Downloading $BIOGRID_URL ..." >&2
    curl -fsSL "$BIOGRID_URL" -o "$BIOGRID_ZIP"
  fi
  echo "Extracting C. elegans tab3 from BioGRID archive ..." >&2
  (cd "$BIOGRID_CACHE" && unzip -o "BIOGRID-ORGANISM-${BIOGRID_VERSION}.tab3.zip" "BIOGRID-ORGANISM-Caenorhabditis_elegans-*.tab3.txt" >/dev/null)
fi

if [[ ! -f "$DB" ]]; then
  echo "Error: $DB not found. Run build_gene_db.R first." >&2
  exit 1
fi

python - "$DB" "$WORMCAT_CSV" "$CACHE/wb.gaf.gz" "$CACHE/go.obo" "$ALLIANCE_ORTHO_TSV_GZ" "$BIOGRID_TAB" "$ALLIANCE_PHENO_JSON_GZ" "$WBPHENO_OBO" "$ANATOMY_OBO" "$ANATOMY_ASSOC" <<'PY'
import csv, gzip, json, sqlite3, sys
from collections import defaultdict
from pathlib import Path

db_path, wormcat_csv, gaf_path, obo_path, ortho_gz, biogrid_tab, pheno_json_gz, wbpheno_obo, anatomy_obo, anatomy_assoc = sys.argv[1:11]
conn = sqlite3.connect(db_path)
conn.executescript("PRAGMA journal_mode=WAL; PRAGMA synchronous=NORMAL;")
cur = conn.cursor()

# ---------- WormCat ----------
cur.executescript("""
DROP TABLE IF EXISTS wormcat;
CREATE TABLE wormcat (
  wormbase_gene TEXT NOT NULL,
  sequence_id   TEXT,
  category_1    TEXT,
  category_2    TEXT,
  category_3    TEXT,
  auto_description TEXT,
  PRIMARY KEY (wormbase_gene)
);
CREATE INDEX idx_wormcat_cat1 ON wormcat(category_1);
CREATE INDEX idx_wormcat_cat2 ON wormcat(category_2);
CREATE INDEX idx_wormcat_cat3 ON wormcat(category_3);
""")

wc_rows = []
if Path(wormcat_csv).is_file():
    with open(wormcat_csv, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            wb = (row.get("Wormbase ID") or "").strip()
            if not wb:
                continue
            def clean(v):
                v = (v or "").strip()
                return None if not v or v.lower() == "unassigned" else v
            wc_rows.append((
                wb,
                (row.get("Sequence ID") or "").strip() or None,
                clean(row.get("Category 1")),
                clean(row.get("Category 2")),
                clean(row.get("Category 3")),
                (row.get("Automated Description") or "").strip() or None,
            ))
    cur.executemany("INSERT OR REPLACE INTO wormcat VALUES (?,?,?,?,?,?)", wc_rows)
    print(f"wormcat: {len(wc_rows)} genes", file=sys.stderr)
else:
    print(f"WARN: wormcat CSV not found at {wormcat_csv}", file=sys.stderr)

# ---------- GO ontology (go.obo) ----------
cur.executescript("""
DROP TABLE IF EXISTS go_terms;
DROP TABLE IF EXISTS go_relations;
DROP TABLE IF EXISTS go_alt_ids;
CREATE TABLE go_terms (
  go_id       TEXT PRIMARY KEY,
  name        TEXT NOT NULL,
  namespace   TEXT NOT NULL,
  definition  TEXT,
  is_obsolete INTEGER NOT NULL DEFAULT 0
);
CREATE TABLE go_relations (
  child_id  TEXT NOT NULL,
  parent_id TEXT NOT NULL,
  rel_type  TEXT NOT NULL,
  PRIMARY KEY (child_id, parent_id, rel_type)
);
CREATE TABLE go_alt_ids (
  alt_id TEXT PRIMARY KEY,
  go_id  TEXT NOT NULL
);
CREATE INDEX idx_gorel_child  ON go_relations(child_id);
CREATE INDEX idx_gorel_parent ON go_relations(parent_id);
""")

NS_MAP = {"biological_process": "BP", "molecular_function": "MF", "cellular_component": "CC"}

terms, rels, alts = [], [], []
cur_term = None

def flush(t):
    if not t or "id" not in t or "name" not in t or "namespace" not in t:
        return
    terms.append((t["id"], t["name"], NS_MAP.get(t["namespace"], t["namespace"]),
                  t.get("def"), 1 if t.get("obsolete") else 0))
    for parent in t.get("is_a", []):
        rels.append((t["id"], parent, "is_a"))
    for rel in t.get("relationship", []):
        # "part_of GO:0008150" → ("part_of", "GO:0008150")
        parts = rel.split(None, 1)
        if len(parts) == 2 and parts[1].startswith("GO:"):
            rels.append((t["id"], parts[1].split()[0], parts[0]))
    for alt in t.get("alt_id", []):
        alts.append((alt, t["id"]))

with open(obo_path, encoding="utf-8") as f:
    in_term = False
    for raw in f:
        line = raw.rstrip("\n")
        if line == "[Term]":
            flush(cur_term)
            cur_term = {"is_a": [], "relationship": [], "alt_id": []}
            in_term = True
            continue
        if line.startswith("["):  # other stanza
            flush(cur_term); cur_term = None; in_term = False
            continue
        if not in_term or not line.strip() or ":" not in line:
            continue
        key, _, val = line.partition(": ")
        val = val.strip()
        if key == "id":
            cur_term["id"] = val
        elif key == "name":
            cur_term["name"] = val
        elif key == "namespace":
            cur_term["namespace"] = val
        elif key == "def":
            # strip surrounding quotes + trailing [refs]
            if val.startswith('"'):
                end = val.rfind('"')
                if end > 0: val = val[1:end]
            cur_term["def"] = val
        elif key == "is_obsolete" and val == "true":
            cur_term["obsolete"] = True
        elif key == "is_a":
            cur_term["is_a"].append(val.split(" !")[0].strip())
        elif key == "relationship":
            cur_term["relationship"].append(val.split(" !")[0].strip())
        elif key == "alt_id":
            cur_term["alt_id"].append(val.strip())
    flush(cur_term)

cur.executemany("INSERT OR IGNORE INTO go_terms VALUES (?,?,?,?,?)", terms)
cur.executemany("INSERT OR IGNORE INTO go_relations VALUES (?,?,?)", rels)
cur.executemany("INSERT OR IGNORE INTO go_alt_ids VALUES (?,?)", alts)
print(f"go_terms: {len(terms)} | go_relations: {len(rels)} | alt_ids: {len(alts)}", file=sys.stderr)

# ---------- GO annotations (wb.gaf.gz) ----------
cur.executescript("""
DROP TABLE IF EXISTS gene_go;
CREATE TABLE gene_go (
  wormbase_gene TEXT NOT NULL,
  go_id         TEXT NOT NULL,
  qualifier     TEXT,
  evidence_code TEXT NOT NULL,
  reference     TEXT,
  with_from     TEXT,
  aspect        TEXT NOT NULL,
  assigned_by   TEXT,
  annotation_date TEXT
);
CREATE INDEX idx_gene_go_gene ON gene_go(wormbase_gene);
CREATE INDEX idx_gene_go_term ON gene_go(go_id);
CREATE INDEX idx_gene_go_evi  ON gene_go(evidence_code);
""")

gaf_rows = []
seen_misses = 0
with gzip.open(gaf_path, "rt", encoding="utf-8") as f:
    for line in f:
        if not line or line.startswith("!"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 15:
            continue
        # GAF 2.2: 0 db, 1 db_object_id, 2 symbol, 3 qualifier, 4 go_id, 5 ref,
        #          6 evidence, 7 with/from, 8 aspect, 9 name, 10 synonyms,
        #          11 type, 12 taxon, 13 date, 14 assigned_by
        db_obj = cols[1].strip()
        if not db_obj.startswith("WBGene"):
            # WormBase GAF uses WBGene IDs in col 1; skip oddities
            seen_misses += 1
            continue
        gaf_rows.append((
            db_obj,
            cols[4].strip(),
            cols[3].strip() or None,
            cols[6].strip(),
            cols[5].strip() or None,
            cols[7].strip() or None,
            cols[8].strip(),
            cols[14].strip() or None,
            cols[13].strip() or None,
        ))

cur.executemany("INSERT INTO gene_go VALUES (?,?,?,?,?,?,?,?,?)", gaf_rows)
print(f"gene_go: {len(gaf_rows)} annotations ({seen_misses} skipped)", file=sys.stderr)

# ---------- Orthologs (Alliance, multi-species) ----------
cur.executescript("""
DROP TABLE IF EXISTS orthologs;
CREATE TABLE orthologs (
  wormbase_gene     TEXT NOT NULL,
  species           TEXT NOT NULL,
  taxon_id          TEXT NOT NULL,
  ortholog_id       TEXT NOT NULL,
  ortholog_symbol   TEXT NOT NULL,
  algorithms        TEXT,
  algorithms_match  INTEGER,
  out_of_algorithms INTEGER,
  is_best_score     INTEGER,
  is_best_rev_score INTEGER,
  PRIMARY KEY (wormbase_gene, taxon_id, ortholog_id)
);
CREATE INDEX idx_ortho_wb       ON orthologs(wormbase_gene);
CREATE INDEX idx_ortho_taxon    ON orthologs(taxon_id);
CREATE INDEX idx_ortho_symbol   ON orthologs(ortholog_symbol);
""")

WORM_TAXON = "NCBITaxon:6239"
ortho_rows = []
seen_keys = set()
if Path(ortho_gz).is_file():
    with gzip.open(ortho_gz, 'rt', encoding='utf-8') as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 13:
                continue
            if cols[0] == 'Gene1ID':  # header line
                continue
            g1_id, g1_sym, g1_tax, g1_sp, g2_id, g2_sym, g2_tax, g2_sp = cols[:8]
            algos = cols[8] or None
            try:
                algos_match = int(cols[9])
            except ValueError:
                algos_match = None
            try:
                out_of = int(cols[10])
            except ValueError:
                out_of = None
            is_best = 1 if cols[11].lower().startswith('y') else 0
            is_best_rev = 1 if cols[12].lower().startswith('y') else 0

            # Keep pairs where one side is C. elegans; normalize so worm is gene1
            if g1_tax == WORM_TAXON:
                wb_id, sp_id, sp_sym, sp_tax, sp_name = g1_id, g2_id, g2_sym, g2_tax, g2_sp
            elif g2_tax == WORM_TAXON:
                wb_id, sp_id, sp_sym, sp_tax, sp_name = g2_id, g1_id, g1_sym, g1_tax, g1_sp
            else:
                continue
            wb_id = wb_id.replace('WB:', '')
            if not wb_id.startswith('WBGene'):
                continue
            key = (wb_id, sp_tax, sp_id)
            if key in seen_keys:
                continue
            seen_keys.add(key)
            ortho_rows.append((wb_id, sp_name, sp_tax, sp_id, sp_sym, algos, algos_match, out_of, is_best, is_best_rev))
    cur.executemany(
        "INSERT OR IGNORE INTO orthologs VALUES (?,?,?,?,?,?,?,?,?,?)",
        ortho_rows
    )
    print(f"orthologs: {len(ortho_rows)} pairs across all species", file=sys.stderr)
else:
    print(f"WARN: ortholog file not found at {ortho_gz}", file=sys.stderr)

# ---------- BioGRID interactions ----------
cur.executescript("""
DROP TABLE IF EXISTS interactions;
CREATE TABLE interactions (
  gene_a              TEXT NOT NULL,
  gene_b              TEXT NOT NULL,
  symbol_a            TEXT,
  symbol_b            TEXT,
  experimental_system TEXT,
  experimental_type   TEXT,
  author              TEXT,
  pmid                TEXT,
  throughput          TEXT,
  score               REAL,
  source              TEXT NOT NULL DEFAULT 'biogrid'
);
CREATE INDEX idx_int_a ON interactions(gene_a);
CREATE INDEX idx_int_b ON interactions(gene_b);
CREATE INDEX idx_int_type ON interactions(experimental_type);
""")

# Build Entrez -> WBGene map from the existing entrez_to_wb table
entrez_to_wb = {}
for entrez, wb in cur.execute("SELECT entrez_id, wormbase_gene FROM entrez_to_wb"):
    if entrez and wb:
        entrez_to_wb[entrez] = wb

ix_rows = []
skipped_no_map = 0
if Path(biogrid_tab).is_file():
    with open(biogrid_tab, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader, None)  # skip header
        for cols in reader:
            if len(cols) < 19:
                continue
            entrez_a = cols[1].strip()
            entrez_b = cols[2].strip()
            wb_a = entrez_to_wb.get(entrez_a)
            wb_b = entrez_to_wb.get(entrez_b)
            if not wb_a or not wb_b:
                skipped_no_map += 1
                continue
            sym_a = cols[7].strip() or None
            sym_b = cols[8].strip() or None
            system = cols[11].strip() or None
            sys_type = cols[12].strip() or None
            author = cols[13].strip() or None
            pmid = cols[14].strip() or None
            throughput = cols[17].strip() or None
            score_raw = cols[18].strip()
            try:
                score = float(score_raw) if score_raw and score_raw != '-' else None
            except ValueError:
                score = None
            ix_rows.append((wb_a, wb_b, sym_a, sym_b, system, sys_type, author, pmid, throughput, score))
    cur.executemany("INSERT INTO interactions VALUES (?,?,?,?,?,?,?,?,?,?,'biogrid')", ix_rows)
    print(f"interactions: {len(ix_rows)} stored ({skipped_no_map} skipped, no Entrez->WB map)", file=sys.stderr)
else:
    print(f"WARN: BioGRID file not found at {biogrid_tab}", file=sys.stderr)

# ---------- Phenotype ontology + gene-phenotype (Alliance + WBPhenotype) ----------
cur.executescript("""
DROP TABLE IF EXISTS phenotype_terms;
DROP TABLE IF EXISTS gene_phenotype;
CREATE TABLE phenotype_terms (
  term_id    TEXT PRIMARY KEY,
  name       TEXT NOT NULL,
  definition TEXT,
  is_obsolete INTEGER NOT NULL DEFAULT 0
);
CREATE TABLE gene_phenotype (
  wormbase_gene TEXT NOT NULL,
  term_id       TEXT NOT NULL,
  statement     TEXT,
  pmid          TEXT,
  source_id     TEXT
);
CREATE INDEX idx_gene_pheno_gene ON gene_phenotype(wormbase_gene);
CREATE INDEX idx_gene_pheno_term ON gene_phenotype(term_id);
""")

# Parse WBPhenotype OBO into phenotype_terms
pheno_terms = []
if Path(wbpheno_obo).is_file():
    cur_term = None
    def flush_p(t):
        if not t or 'id' not in t or 'name' not in t:
            return
        if not t['id'].startswith('WBPhenotype:'):
            return
        pheno_terms.append((t['id'], t['name'], t.get('def'), 1 if t.get('obsolete') else 0))
    with open(wbpheno_obo, encoding='utf-8') as f:
        in_term = False
        for raw in f:
            line = raw.rstrip('\n')
            if line == '[Term]':
                flush_p(cur_term)
                cur_term = {}
                in_term = True
                continue
            if line.startswith('['):
                flush_p(cur_term); cur_term = None; in_term = False
                continue
            if not in_term or not line.strip() or ':' not in line:
                continue
            k, _, v = line.partition(': ')
            v = v.strip()
            if k == 'id': cur_term['id'] = v
            elif k == 'name': cur_term['name'] = v
            elif k == 'def':
                if v.startswith('"'):
                    end = v.rfind('"')
                    if end > 0: v = v[1:end]
                cur_term['def'] = v
            elif k == 'is_obsolete' and v == 'true':
                cur_term['obsolete'] = True
        flush_p(cur_term)
    cur.executemany("INSERT OR IGNORE INTO phenotype_terms VALUES (?,?,?,?)", pheno_terms)
    print(f"phenotype_terms: {len(pheno_terms)}", file=sys.stderr)

# Parse Alliance PHENOTYPE_WB.json.gz, gene-level entries only
gp_rows = []
seen_gp = set()
if Path(pheno_json_gz).is_file():
    with gzip.open(pheno_json_gz, 'rt', encoding='utf-8') as f:
        payload = json.load(f)
    for e in payload.get('data', []):
        oid = e.get('objectId', '')
        if not oid.startswith('WB:WBGene'):
            continue
        wb = oid.replace('WB:', '')
        statement = (e.get('phenotypeStatement') or '').strip() or None
        pmid = (e.get('evidence') or {}).get('publicationId')
        for term in e.get('phenotypeTermIdentifiers', []):
            tid_raw = (term.get('termId') or '').strip()
            if not tid_raw:
                continue
            # Normalize "WB:WBPhenotype:0000505" -> "WBPhenotype:0000505"
            tid = tid_raw.replace('WB:', '')
            if not tid.startswith('WBPhenotype:'):
                continue
            key = (wb, tid, statement, pmid)
            if key in seen_gp:
                continue
            seen_gp.add(key)
            gp_rows.append((wb, tid, statement, pmid, oid))
    cur.executemany("INSERT INTO gene_phenotype VALUES (?,?,?,?,?)", gp_rows)
    print(f"gene_phenotype: {len(gp_rows)} gene-level associations", file=sys.stderr)

# ---------- Anatomy ontology + gene-anatomy (TEA, WormBase) ----------
cur.executescript("""
DROP TABLE IF EXISTS anatomy_terms;
DROP TABLE IF EXISTS gene_anatomy;
CREATE TABLE anatomy_terms (
  term_id    TEXT PRIMARY KEY,
  name       TEXT NOT NULL,
  definition TEXT,
  is_obsolete INTEGER NOT NULL DEFAULT 0
);
CREATE TABLE gene_anatomy (
  wormbase_gene TEXT NOT NULL,
  term_id       TEXT NOT NULL,
  evidence_code TEXT,
  PRIMARY KEY (wormbase_gene, term_id, evidence_code)
);
CREATE INDEX idx_gene_anat_gene ON gene_anatomy(wormbase_gene);
CREATE INDEX idx_gene_anat_term ON gene_anatomy(term_id);
""")

# Parse WBbt OBO into anatomy_terms (same logic as WBPhenotype OBO above).
anat_terms = []
if Path(anatomy_obo).is_file():
    cur_term = None
    def flush_a(t):
        if not t or 'id' not in t or 'name' not in t:
            return
        if not t['id'].startswith('WBbt:'):
            return
        anat_terms.append((t['id'], t['name'], t.get('def'), 1 if t.get('obsolete') else 0))
    with open(anatomy_obo, encoding='utf-8') as f:
        in_term = False
        for raw in f:
            line = raw.rstrip('\n')
            if line == '[Term]':
                flush_a(cur_term)
                cur_term = {}
                in_term = True
                continue
            if line.startswith('['):
                flush_a(cur_term); cur_term = None; in_term = False
                continue
            if not in_term or not line.strip() or ':' not in line:
                continue
            k, _, v = line.partition(': ')
            v = v.strip()
            if k == 'id': cur_term['id'] = v
            elif k == 'name': cur_term['name'] = v
            elif k == 'def':
                if v.startswith('"'):
                    end = v.rfind('"')
                    if end > 0: v = v[1:end]
                cur_term['def'] = v
            elif k == 'is_obsolete' and v == 'true':
                cur_term['obsolete'] = True
        flush_a(cur_term)
    cur.executemany("INSERT OR IGNORE INTO anatomy_terms VALUES (?,?,?,?)", anat_terms)
    print(f"anatomy_terms: {len(anat_terms)}", file=sys.stderr)
else:
    print(f"WARN: anatomy OBO not found at {anatomy_obo}", file=sys.stderr)

# Parse Alliance EXPRESSION-ALLIANCE_WB.tsv. Header line starts with "Species".
# Columns (0-indexed):
#   2 = GeneID ("WB:WBGene...")
#   7 = AssayTermName (we keep as evidence_code so users can filter)
#  16 = AnatomyTermID ("WBbt:...")
#  18 = AnatomyTermQualifierIDs ("not" anywhere => drop)
ga_rows = []
seen_ga = set()
if Path(anatomy_assoc).is_file():
    open_fn = gzip.open if anatomy_assoc.endswith('.gz') else open
    with open_fn(anatomy_assoc, 'rt', encoding='utf-8', errors='replace') as f:
        in_header_block = True
        for raw in f:
            if not raw or raw.startswith('#'):
                continue
            parts = raw.rstrip('\n').split('\t')
            if in_header_block:
                if parts and parts[0] == 'Species':
                    in_header_block = False
                continue
            if len(parts) < 18:
                continue
            wb_raw = parts[2].strip()
            assay = parts[7].strip() if len(parts) > 7 else ''
            term = parts[16].strip()
            qualifier = parts[18].strip().lower() if len(parts) > 18 else ''
            if not wb_raw.startswith('WB:WBGene'):
                continue
            wb = wb_raw[3:]  # strip "WB:"
            if not term.startswith('WBbt:'):
                continue
            if 'not' in qualifier.split('|'):
                continue
            key = (wb, term, assay)
            if key in seen_ga:
                continue
            seen_ga.add(key)
            ga_rows.append((wb, term, assay or None))
    cur.executemany("INSERT OR IGNORE INTO gene_anatomy VALUES (?,?,?)", ga_rows)
    print(f"gene_anatomy: {len(ga_rows)} associations", file=sys.stderr)
else:
    print(f"WARN: anatomy associations not found at {anatomy_assoc}", file=sys.stderr)

# ---------- GO IDF (per-term inverse document frequency) ----------
cur.executescript("""
DROP TABLE IF EXISTS go_term_freq;
CREATE TABLE go_term_freq (
  go_id TEXT PRIMARY KEY,
  df    INTEGER NOT NULL,
  idf   REAL NOT NULL
);
""")

import math
N_total = cur.execute("SELECT COUNT(DISTINCT wormbase_gene) FROM gene_go").fetchone()[0] or 1
idf_rows = []
for go_id, df in cur.execute("SELECT go_id, COUNT(DISTINCT wormbase_gene) AS df FROM gene_go GROUP BY go_id"):
    if df > 0:
        idf = math.log(N_total / df)
        idf_rows.append((go_id, df, max(idf, 0.0)))
cur.executemany("INSERT INTO go_term_freq VALUES (?,?,?)", idf_rows)
print(f"go_term_freq: {len(idf_rows)} terms (N_total={N_total})", file=sys.stderr)

# ---------- FTS5 indexes ----------
cur.executescript("""
DROP TABLE IF EXISTS genes_fts;
CREATE VIRTUAL TABLE genes_fts USING fts5(
  wormbase_gene UNINDEXED,
  gene_name, sequence_name, other_names,
  gene_class_description, concise_description, legacy_description,
  tokenize = "unicode61 remove_diacritics 2"
);
INSERT INTO genes_fts(wormbase_gene, gene_name, sequence_name, other_names,
                     gene_class_description, concise_description, legacy_description)
SELECT wormbase_gene,
       COALESCE(gene_name,''),
       COALESCE(sequence_name,''),
       COALESCE(other_names,''),
       COALESCE(gene_class_description,''),
       COALESCE(concise_description,''),
       COALESCE(legacy_description,'')
FROM annotations;

DROP TABLE IF EXISTS go_fts;
CREATE VIRTUAL TABLE go_fts USING fts5(
  go_id UNINDEXED, name, definition, namespace UNINDEXED,
  tokenize = "unicode61 remove_diacritics 2"
);
INSERT INTO go_fts(go_id, name, definition, namespace)
SELECT go_id, name, COALESCE(definition,''), namespace
FROM go_terms WHERE is_obsolete = 0;

DROP TABLE IF EXISTS wormcat_fts;
CREATE VIRTUAL TABLE wormcat_fts USING fts5(
  wormbase_gene UNINDEXED, category_1, category_2, category_3, auto_description,
  tokenize = "unicode61 remove_diacritics 2"
);
INSERT INTO wormcat_fts(wormbase_gene, category_1, category_2, category_3, auto_description)
SELECT wormbase_gene, COALESCE(category_1,''), COALESCE(category_2,''),
       COALESCE(category_3,''), COALESCE(auto_description,'')
FROM wormcat;
""")
print("fts5: genes_fts, go_fts, wormcat_fts built", file=sys.stderr)

conn.commit()
conn.executescript("VACUUM; ANALYZE;")
conn.close()
print(f"OK -> {db_path}", file=sys.stderr)
PY

echo
echo "Tables in $DB:"
sqlite3 "$DB" ".tables"
echo
echo "Sizes:"
sqlite3 "$DB" "SELECT 'annotations', COUNT(*) FROM annotations
UNION ALL SELECT 'wormcat', COUNT(*) FROM wormcat
UNION ALL SELECT 'go_terms', COUNT(*) FROM go_terms
UNION ALL SELECT 'go_relations', COUNT(*) FROM go_relations
UNION ALL SELECT 'gene_go', COUNT(*) FROM gene_go
UNION ALL SELECT 'orthologs', COUNT(*) FROM orthologs
UNION ALL SELECT 'interactions', COUNT(*) FROM interactions
UNION ALL SELECT 'phenotype_terms', COUNT(*) FROM phenotype_terms
UNION ALL SELECT 'gene_phenotype', COUNT(*) FROM gene_phenotype;"
echo
echo "DB size: $(du -h "$DB" | cut -f1)"
