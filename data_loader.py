"""Data access layer for the poplar / Arabidopsis activation-domain viewer.

Coordinate convention used everywhere in this module: 1-based, inclusive
(start, end), matching Pfam/HMMER env coordinates and the tile `position`
column in the prediction pickles.
"""
from __future__ import annotations

import csv
import io
import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

DATA_DIR = Path(__file__).parent / "data_for_viewer"
AD_THRESHOLD = 1.0
AD_MIN_LEN = 10
TILE_LEN = 40

_PRED_COLS = ["family", "protein_ID", "gene_ID", "seq", "len", "preds",
              "disorder", "position", "AAseq", "activity", "activity_avg"]

# Per-species TF tables with CDS (first existing file wins; looked up in the data
# directory, then the project root). `family` = whether its Family column overrides.
TF_TABLES = {
    "Poplar": {"files": ["Poplar_filteredTFsClean.xlsx", "Poplar_filteredTFsClean.csv"], "family": True},
    "Arabidopsis": {"files": ["AthAllTFCDS.xlsx", "AthAllTFCDS.csv"], "family": False},
}
BULK_SUMMARY_COLS = ["Gene_ID", "Family", "Has_AD", "Num_ADs", "Peak_Predicted_Score", "Species",
                     "AD_Coordinates", "AD_Peak_Scores", "Pfam_Domains", "Empirical_Active_Tiles",
                     "Protein_ID", "Length", "CDS_Available"]
BULK_SEQUENCE_COLS = ["Full_Protein_Seq", "Full_CDS_Seq", "AD_Peptide_Seqs", "AD_CDS_Seqs",
                      "Pfam_Peptide_Seqs", "Pfam_CDS_Seqs"]
FASTA_KINDS = {"ad_aa": "AD peptides", "ad_cds": "AD CDS", "protein": "full proteins", "cds": "full CDS"}

_ID_COLS = ["locus", "gene_id", "agi", "protein_id", "gene", "id"]
_CDS_COLS = ["cds", "sequence", "cds_sequence", "seq", "nt_seq"]
EMPIRICAL_SOURCE = {"Poplar": "This study", "Arabidopsis": "Morffy et al., 2024"}
MORFFY_URL = "https://doi.org/10.1038/s41586-024-07707-3"
MORFFY_CITATION = ("Morffy, N. et al. Identification of plant transcriptional activation domains. "
                   "*Nature* **632**, 166–173 (2024).")

# Pfam models that describe the same domain (no clan table ships with the data, so
# this is a curated stand-in). Overlapping/touching hits within one group are merged.
PFAM_GROUPS = {
    "MYB DBD": ["Myb_DNA-binding", "Myb_DNA-bind_3", "Myb_DNA-bind_4", "Myb_DNA-bind_5",
                "Myb_DNA-bind_6", "Myb_DNA-bind_7", "MADF_DNA_bdg", "SANT"],
    "bZIP DBD": ["bZIP_1", "bZIP_2", "bZIP_Maf"],
    "Homeodomain": ["Homeodomain", "Homeobox_KN", "Homez"],
    "C2H2 zinc finger": ["zf-C2H2", "zf-C2H2_4", "zf-C2H2_6", "zf-C2H2_jaz", "zf-H2C2_2", "zf-met"],
    "CCCH zinc finger": ["zf-CCCH", "zf-CCCH_2", "zf-CCCH_3", "zf-CCCH_4", "zf_CCCH_4"],
    "RING finger": ["zf-RING_2", "zf-C3HC4", "zf-C3HC4_2", "zf-C3HC4_3", "zf-RING_UBOX"],
    "Ankyrin repeats": ["Ank", "Ank_2", "Ank_3", "Ank_4", "Ank_5"],
    "Histone fold (NF-Y)": ["Histone", "CBFD_NFYB_HMF"],
}
_PFAM_GROUP_OF = {m: g for g, ms in PFAM_GROUPS.items() for m in ms}


@dataclass
class TFRecord:
    species: str            # "Poplar" | "Arabidopsis"
    gene_id: str
    protein_id: str
    family: str
    seq: str
    preds: np.ndarray       # per-residue predicted activity, len == len(seq)
    disorder: np.ndarray    # per-residue disorder score, len == len(seq)
    activity_avg: np.ndarray  # per-residue mean of overlapping empirical tiles
    tiles: pd.DataFrame     # start, end, activity, tile_seq, placement
    pfam: pd.DataFrame      # merged domains: pfam_name (label), pfam_acc, start, end, i_evalue, is_dbd, members
    pfam_raw: pd.DataFrame  # unmerged hmmscan hits
    cds: str | None = None  # coding sequence (incl. stop codon if present)
    ads: pd.DataFrame = field(default_factory=pd.DataFrame)

    @property
    def length(self) -> int:
        return len(self.seq)

    @property
    def empirical_source(self) -> str:
        return EMPIRICAL_SOURCE[self.species]

    def subseq(self, start, end) -> str:
        """1-based inclusive slice of the protein sequence."""
        start_idx, end_idx = _coords(start, end)
        return self.seq[start_idx - 1:end_idx]

    def cds_slice(self, start, end) -> str | None:
        """Nucleotides encoding residues start..end (1-based, inclusive):
        cds[(start - 1) * 3 : end * 3]. None if no CDS is available."""
        if self.cds is None:
            return None
        start_idx, end_idx = _coords(start, end)
        start_nt = (start_idx - 1) * 3
        end_nt = end_idx * 3
        return self.cds[start_nt:end_nt]


class TFData:
    """Loads all files once; `get(gene_or_protein_id)` returns a TFRecord."""

    def __init__(self, data_dir: Path | str = DATA_DIR):
        d = Path(data_dir)
        with warnings.catch_warnings():  # numpy.core pickle deprecation noise
            warnings.simplefilter("ignore", DeprecationWarning)
            ptr = pd.read_pickle(d / "ptr_all_preds.pkl")[_PRED_COLS]
            ath = pd.read_pickle(d / "ath_all_preds.pkl")[_PRED_COLS]
            self.mapping = pd.read_pickle(d / "mapping.pkl")
        ptr["species"] = "Poplar"
        ath["species"] = "Arabidopsis"
        self.preds = pd.concat([ptr, ath], ignore_index=True)

        self._records: dict[int, TFRecord] = {}
        # Case-insensitive lookup by gene ID or protein ID -> row index
        self._lookup: dict[str, int] = {}
        for i, (g, p) in enumerate(zip(self.preds.gene_ID, self.preds.protein_ID)):
            self._lookup[g.upper()] = i
            self._lookup[p.upper()] = i

        pf = pd.read_csv(d / "pfam_hits_all_proteins.csv")
        pf = pf.rename(columns={"env_start": "start", "env_end": "end"})
        self._pfam = {k: g.drop(columns="protein_id").sort_values("start").reset_index(drop=True)
                      for k, g in pf.groupby("protein_id")}

        self._load_tf_table(d)
        self._build_orthogroups(d)

    # ---- per-species TF tables (CDS, optionally family) -----------------------
    def _load_tf_table(self, d: Path) -> None:
        """Load CDS (and poplar Family) from the TF_TABLES files.

        Rows are matched to proteins by protein ID (e.g. AT4G29930.1); when the
        exact isoform is absent, any isoform of the same gene whose CDS
        translates to the stored protein is used. A CDS is kept only if it
        encodes exactly the stored protein (an optional trailing stop codon is
        allowed), so AA->nt slicing is always valid.
        """
        self._cds: dict[str, str] = {}
        self.cds_rejected: list[str] = []
        self.cds_sources: dict[str, Path] = {}
        self.load_notes: list[str] = []
        fam = self.preds.family.copy()
        for species, spec in TF_TABLES.items():
            tbl, path = _read_tf_table(spec["files"], [d, d.parent], self.load_notes)
            if tbl is None:
                continue
            self.cds_sources[species] = path
            by_id = tbl.set_index("_id")
            isoforms: dict[str, list[str]] = {}
            for loc in tbl["_id"]:
                isoforms.setdefault(loc.upper().rsplit(".", 1)[0], []).append(loc)
            for i in np.flatnonzero((self.preds.species == species).values):
                pid, gid = self.preds.protein_ID.iat[i], self.preds.gene_ID.iat[i]
                seq = self.preds.seq.iat[i]
                cands = [pid] if pid in by_id.index else isoforms.get(gid.upper(), [])
                hit = next((c for c in cands if _encodes(by_id.at[c, "_cds"], seq)), None)
                if hit is None:
                    if cands:
                        self.cds_rejected.append(pid)
                    continue
                self._cds[pid] = by_id.at[hit, "_cds"]
                f = by_id.at[hit, "_family"] if "_family" in by_id.columns else None
                if spec["family"] and isinstance(f, str) and f:
                    fam.iat[i] = f
        self.preds["family"] = fam

    def cds_coverage(self) -> dict[str, tuple[int, int]]:
        """species -> (TFs with CDS, total TFs)."""
        has = self.preds.protein_ID.isin(self._cds.keys())
        return {sp: (int(has[self.preds.species == sp].sum()), int((self.preds.species == sp).sum()))
                for sp in ["Poplar", "Arabidopsis"]}

    # ---- families ----------------------------------------------------------
    def family_counts(self, species: str | None = None) -> dict[str, int]:
        df = self.preds if species is None else self.preds[self.preds.species == species]
        return dict(sorted(df.family.value_counts().items(), key=lambda kv: kv[0].lower()))

    # ---- orthogroups -------------------------------------------------------
    def _build_orthogroups(self, d: Path) -> None:
        """Orthogroups = connected components of mapping.pkl restricted to TF-TF edges.

        This reproduces exactly the 641 groups that have precomputed MSAs in
        ad_locus_conservation_raw.pkl (the remaining components lack an MSA and
        are shown on residue coordinates).
        """
        m = self.mapping.rename(columns={"Gene Name": "ptr", "Ortholog Gene Name": "ath",
                                         "Relationship": "relationship"})
        tf = set(self.preds.gene_ID)
        is_tf_edge = m.ptr.isin(tf) & m.ath.isin(tf)
        self.tf_edges = m[is_tf_edge].reset_index(drop=True)
        # edges from a TF to a gene outside the TF set (no prediction data)
        self.external_edges = m[~is_tf_edge & (m.ptr.isin(tf) | m.ath.isin(tf))].reset_index(drop=True)

        parent: dict[str, str] = {}

        def find(x):
            parent.setdefault(x, x)
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for a, b in zip(self.tf_edges.ptr, self.tf_edges.ath):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
        comps: dict[str, list[str]] = {}
        for g in parent:
            comps.setdefault(find(g), []).append(g)
        self._group_of: dict[str, frozenset] = {}
        for members in comps.values():
            fs = frozenset(members)
            for g in members:
                self._group_of[g] = fs

        # Precomputed MSAs, keyed by member gene set -> {gene_id: aligned_seq}
        self._msa: dict[frozenset, dict[str, str]] = {}
        path = d / "ad_locus_conservation_raw.pkl"
        if path.exists():
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", DeprecationWarning)
                raw = pd.read_pickle(path)
            gene_of = dict(zip(self.preds.protein_ID, self.preds.gene_ID))
            for grp in raw:
                aln = {gene_of.get(pid, pid.rsplit(".", 1)[0]): v["aligned_seq"] for pid, v in grp.items()}
                self._msa[frozenset(aln)] = aln

    def orthogroup(self, gene_or_protein_id: str) -> "Orthogroup":
        gid = self.get(gene_or_protein_id).gene_id
        members = self._group_of.get(gid, frozenset([gid]))
        species = dict(zip(self.preds.gene_ID, self.preds.species))
        ptr = sorted(g for g in members if species[g] == "Poplar")
        ath = sorted(g for g in members if species[g] == "Arabidopsis")
        edges = self.tf_edges[self.tf_edges.ptr.isin(members)].reset_index(drop=True)
        ext = self.external_edges
        ext = ext[ext.ptr.isin(members) | ext.ath.isin(members)].reset_index(drop=True)
        msa = self._msa.get(members)
        if msa is not None:  # sanity: alignment must degap to the stored sequence
            for g, a in msa.items():
                if a.replace("-", "") != self.get(g).seq:
                    msa = None
                    break
        return Orthogroup(query=gid, ptr=ptr, ath=ath, edges=edges, external_edges=ext, msa=msa)

    # ---- bulk lookup -------------------------------------------------------
    def resolve(self, tokens: list[str]) -> tuple[list[str], list[str]]:
        """Map raw tokens to unique gene IDs (input order kept).

        Returns (gene_ids, not_found). Gene and protein IDs for the same gene
        collapse to one entry; matching is case-insensitive.
        """
        found, missing, seen, seen_missing = [], [], set(), set()
        for t in tokens:
            key = t.strip().strip(_EDGE_PUNCT).upper()
            if not key:
                continue
            if key not in self._lookup:  # other transcript, e.g. AT1G71930.2 -> AT1G71930
                key = re.sub(r"\.\d+$", "", key)
            if key in self._lookup:
                gid = self.preds.gene_ID.iat[self._lookup[key]]
                if gid not in seen:
                    seen.add(gid)
                    found.append(gid)
            elif key not in seen_missing:
                seen_missing.add(key)
                missing.append(t.strip())
        return found, missing

    def bulk_summary(self, gene_ids: list[str]) -> pd.DataFrame:
        """One row per gene: summary columns first, then full sequence / slice columns."""
        rows = []
        for g in gene_ids:
            rec = self.get(g)
            ads = rec.ads
            ad_iv = [(*_coords(s, e), f"AD{i}") for i, (s, e) in enumerate(zip(ads.Start, ads.End), 1)]
            pf_iv = [(*_coords(s, e), n) for s, e, n in zip(rec.pfam.start, rec.pfam.end, rec.pfam.pfam_name)]
            n_tiles_active = int((rec.tiles.activity >= AD_THRESHOLD).sum()) if len(rec.tiles) else None
            rows.append({
                "Gene_ID": rec.gene_id,
                "Family": rec.family,
                "Has_AD": bool(len(ads)),
                "Num_ADs": int(len(ads)),
                "Peak_Predicted_Score": round(float(np.max(rec.preds)), 3),
                "Species": "Ptr" if rec.species == "Poplar" else "Ath",
                "AD_Coordinates": "; ".join(f"{s}-{e}" for s, e, _ in ad_iv) or "None",
                "AD_Peak_Scores": "; ".join(f"{p:.2f}" for p in ads["Peak Score"]) or "None",
                "Pfam_Domains": ", ".join(rec.pfam.pfam_name) or "None",
                "Empirical_Active_Tiles": n_tiles_active if n_tiles_active is not None else "N/A",
                "Protein_ID": rec.protein_id,
                "Length": rec.length,
                "CDS_Available": rec.cds is not None,
                "Full_Protein_Seq": rec.seq,
                "Full_CDS_Seq": rec.cds if rec.cds is not None else "N/A",
                "AD_Peptide_Seqs": _labelled(ad_iv, rec.subseq),
                "AD_CDS_Seqs": _labelled(ad_iv, rec.cds_slice),
                "Pfam_Peptide_Seqs": _labelled(pf_iv, rec.subseq),
                "Pfam_CDS_Seqs": _labelled(pf_iv, rec.cds_slice),
            })
        return pd.DataFrame(rows, columns=BULK_SUMMARY_COLS + BULK_SEQUENCE_COLS)

    def fasta(self, gene_ids: list[str], *, what: str, width: int = 60) -> str:
        """FASTA export for the queried genes (AA coordinates, 1-based inclusive).

        what="ad_aa"   >GeneID_AD1_start-end peak_score=X.XX      (AD peptides)
        what="ad_cds"  >GeneID_AD1_CDS_start-end peak_score=X.XX  (cds[(start-1)*3 : end*3])
        what="protein" >GeneID protein_id species family length=N aa
        what="cds"     >GeneID protein_id species family length=N nt
        Records whose CDS is unavailable are omitted from the CDS exports.
        Arguments after gene_ids are keyword-only.
        """
        if what not in FASTA_KINDS:
            raise ValueError(f"what must be one of {list(FASTA_KINDS)}, got {what!r}")
        width = int(width)
        out = []

        def emit(header, seq):
            if seq:
                out.append(">" + header)
                out.extend(seq[i:i + width] for i in range(0, len(seq), width))

        for g in gene_ids:
            rec = self.get(g)
            meta = f"{rec.protein_id} {rec.species} {rec.family}"
            if what == "protein":
                emit(f"{rec.gene_id} {meta} length={rec.length} aa", rec.seq)
            elif what == "cds":
                if rec.cds is not None:
                    emit(f"{rec.gene_id} {meta} length={len(rec.cds)} nt", rec.cds)
            else:
                for i, (s, e, pk) in enumerate(zip(rec.ads.Start, rec.ads.End, rec.ads["Peak Score"]), 1):
                    s, e = _coords(s, e)
                    if what == "ad_aa":
                        emit(f"{rec.gene_id}_AD{i}_{s}-{e} peak_score={float(pk):.2f}", rec.subseq(s, e))
                    else:
                        emit(f"{rec.gene_id}_AD{i}_CDS_{s}-{e} peak_score={float(pk):.2f}", rec.cds_slice(s, e))
        return "\n".join(out) + ("\n" if out else "")

    # ---- full-dataset tables (About-page downloads) --------------------------
    def orthology_table(self) -> pd.DataFrame:
        """Every mapping.pkl row, annotated with TF status and its TF orthogroup.

        Orthogroup IDs (OG0001, ...) are assigned deterministically to the TF-TF
        connected components used by the Orthogroup Viewer, ordered by their
        sorted member list; rows involving a non-TF gene have no orthogroup.
        """
        groups = sorted({tuple(sorted(g)) for g in self._group_of.values()})
        og_id = {frozenset(g): f"OG{i:04d}" for i, g in enumerate(groups, 1)}
        species = dict(zip(self.preds.gene_ID, self.preds.species))
        info = {}
        for g in groups:
            fs = frozenset(g)
            n_ptr = sum(species[x] == "Poplar" for x in g)
            info[fs] = (og_id[fs], n_ptr, len(g) - n_ptr, fs in self._msa)
        m = self.mapping.rename(columns={"Gene Name": "Poplar_Gene", "Ortholog Gene Name": "Arabidopsis_Gene",
                                         "Relationship": "Relationship"}).reset_index(drop=True)
        tf = set(self.preds.gene_ID)
        m["Poplar_is_TF"] = m.Poplar_Gene.isin(tf)
        m["Arabidopsis_is_TF"] = m.Arabidopsis_Gene.isin(tf)
        both = m.Poplar_is_TF & m.Arabidopsis_is_TF
        rows = [info[self._group_of[p]] if b else ("", None, None, None)
                for p, b in zip(m.Poplar_Gene, both)]
        m[["Orthogroup_ID", "Orthogroup_Ptr_Genes", "Orthogroup_Ath_Genes", "MSA_Available"]] = pd.DataFrame(
            rows, index=m.index)
        m["Copy_Ratio_Ptr_Ath"] = [f"{int(p)}:{int(a)}" if b else "" for p, a, b in
                                   zip(m.Orthogroup_Ptr_Genes.fillna(0), m.Orthogroup_Ath_Genes.fillna(0), both)]
        for c in ["Orthogroup_Ptr_Genes", "Orthogroup_Ath_Genes"]:
            m[c] = m[c].astype("Int64")
        return m.sort_values(["Orthogroup_ID", "Poplar_Gene", "Arabidopsis_Gene"],
                             key=lambda s: s.replace("", "~") if s.name == "Orthogroup_ID" else s,
                             ).reset_index(drop=True)

    def pfam_table(self) -> pd.DataFrame:
        """Merged Pfam domain intervals for every TF (same merging as the viewers)."""
        frames = []
        for pid, gid, sp in zip(self.preds.protein_ID, self.preds.gene_ID, self.preds.species):
            if pid not in self._pfam:
                continue
            d = self.get(pid).pfam.copy()  # memoized record: shares work with the matrix export
            d.insert(0, "Species", "Ptr" if sp == "Poplar" else "Ath")
            d.insert(0, "Protein_ID", pid)
            d.insert(0, "Gene_ID", gid)
            frames.append(d)
        out = pd.concat(frames, ignore_index=True).rename(columns={
            "pfam_name": "Domain", "start": "Start", "end": "End", "is_dbd": "DNA_Binding",
            "members": "Pfam_Models", "pfam_acc": "Pfam_Accessions", "i_evalue": "Best_i_Evalue",
            "n_hits": "N_Hits"})
        out.insert(out.columns.get_loc("End") + 1, "Length", out.End - out.Start + 1)
        out["DNA_Binding"] = out.DNA_Binding.astype(bool)
        return out[["Gene_ID", "Protein_ID", "Species", "Domain", "Start", "End", "Length", "DNA_Binding",
                    "Pfam_Models", "Pfam_Accessions", "Best_i_Evalue", "N_Hits"]]

    def orthogroup_summary(self, og: "Orthogroup") -> pd.DataFrame:
        rows = []
        for g in og.members:
            rec = self.get(g)
            ads = rec.ads
            rows.append({
                "Gene ID": g,
                "Species": rec.species,
                "Family": rec.family,
                "Length (aa)": rec.length,
                "Activator": "Yes" if len(ads) else "No",
                "# ADs": len(ads),
                "AD coordinates": "; ".join(f"{s}–{e}" for s, e in zip(ads.Start, ads.End)),
                "AD peak scores": "; ".join(f"{p:.2f}" for p in ads["Peak Score"]),
                "Max pred": round(float(rec.preds.max()), 2),
            })
        return pd.DataFrame(rows)

    # ---- listing -----------------------------------------------------------
    def gene_ids(self, species: str | None = None, family: str | None = None) -> list[str]:
        df = self.preds
        if species is not None:
            df = df[df.species == species]
        if family is not None:
            df = df[df.family == family]
        return sorted(df.gene_ID.tolist())

    def has(self, gene_or_protein_id: str) -> bool:
        return gene_or_protein_id.strip().upper() in self._lookup

    # ---- single record -----------------------------------------------------
    def get(self, gene_or_protein_id: str) -> TFRecord:
        key = gene_or_protein_id.strip().upper()
        if key not in self._lookup:
            raise KeyError(f"Unknown gene/protein ID: {gene_or_protein_id!r}")
        idx = self._lookup[key]
        if idx in self._records:
            return self._records[idx]
        r = self.preds.iloc[idx]
        seq = r.seq
        preds = np.asarray(r.preds, dtype=float)
        rec = TFRecord(
            species=r.species,
            gene_id=r.gene_ID,
            protein_id=r.protein_ID,
            family=r.family,
            seq=seq,
            preds=preds,
            disorder=np.asarray(r.disorder, dtype=float),
            activity_avg=np.asarray(r.activity_avg, dtype=float),
            tiles=_place_tiles(seq, r.position, r.AAseq, r.activity),
            pfam=merge_pfam(self._pfam.get(r.protein_ID, _empty_pfam())),
            pfam_raw=self._pfam.get(r.protein_ID, _empty_pfam()),
            cds=self._cds.get(r.protein_ID),
        )
        rec.ads = call_ads(preds, seq)
        self._records[idx] = rec  # memoize; records are treated as read-only by callers
        return rec


@dataclass
class Orthogroup:
    query: str
    ptr: list[str]                 # poplar gene IDs in the group
    ath: list[str]                 # Arabidopsis gene IDs in the group
    edges: pd.DataFrame            # TF-TF mapping rows within the group (ptr, ath, relationship)
    external_edges: pd.DataFrame   # mapping rows to genes outside the TF dataset
    msa: dict[str, str] | None     # gene_id -> aligned sequence, if a precomputed MSA exists

    @property
    def members(self) -> list[str]:
        return self.ptr + self.ath

    @property
    def ratio(self) -> str:
        return f"{len(self.ptr)}:{len(self.ath)}"

    @property
    def aligned(self) -> bool:
        return self.msa is not None

    @property
    def aln_len(self) -> int | None:
        return len(next(iter(self.msa.values()))) if self.msa else None

    def columns(self, gene_id: str) -> np.ndarray | None:
        """1-based alignment column for each residue of `gene_id` (None if no MSA)."""
        if not self.msa:
            return None
        a = np.frombuffer(self.msa[gene_id].encode(), dtype=np.uint8)
        return np.flatnonzero(a != ord("-")) + 1


# ---- helpers ---------------------------------------------------------------
_SPLIT = re.compile(r"[\s,;|]+")   # \s covers newlines, tabs and unicode spaces
_EDGE_PUNCT = "\"'`()[]{}<>.,;:!?*"


def parse_id_text(text: str) -> list[str]:
    """Split pasted text on any mix of newlines, commas, tabs, semicolons, pipes and
    spaces; strip quotes, brackets and trailing punctuation from each token."""
    toks = (t.strip(_EDGE_PUNCT) for t in _SPLIT.split(text or ""))
    return [t for t in toks if t]


_ID_LIKE = re.compile(r"^\s*[\"']?(Potri\.|AT[1-5CM]G\d)", re.IGNORECASE)


def parse_id_file(raw: bytes, filename: str = "") -> list[str]:
    """Tokens from an uploaded .txt/.csv.

    .csv: uses the column with the most gene-ID-like cells (so headers and
    annotation columns are ignored). Anything else: treated as free text.
    """
    text = raw.decode("utf-8-sig", errors="replace")
    if filename.lower().endswith(".csv"):
        rows = [r for r in csv.reader(io.StringIO(text)) if any(c.strip() for c in r)]
        ncol = max((len(r) for r in rows), default=0)
        if ncol > 1:
            hits = [sum(bool(_ID_LIKE.match(r[j])) for r in rows if j < len(r)) for j in range(ncol)]
            j = int(np.argmax(hits))
            if hits[j]:
                return [r[j].strip().strip("\"'") for r in rows
                        if j < len(r) and _ID_LIKE.match(r[j])]
        # single column: drop a header row that doesn't look like an ID
        if rows and not _ID_LIKE.match(rows[0][0]):
            rows = rows[1:]
        return parse_id_text("\n".join(",".join(r) for r in rows))
    return parse_id_text(text)


def _labelled(intervals, slicer) -> str:
    """'AD1 (120-155): SEQ; AD2 (...): SEQ' for (start, end, label) intervals; 'None' if empty,
    'N/A' if the slicer has no sequence (no CDS)."""
    if not intervals:
        return "None"
    parts = []
    for s, e, label in intervals:
        seq = slicer(s, e)
        if seq is None:
            return "N/A"
        parts.append(f"{label} ({s}-{e}): {seq}")
    return "; ".join(parts)


def _coords(start, end) -> tuple[int, int]:
    """Cast 1-based inclusive coordinates to int (accepts numpy ints, floats, numeric strings)."""
    start_idx, end_idx = int(float(start)), int(float(end))
    if start_idx < 1 or end_idx < start_idx:
        raise ValueError(f"invalid 1-based interval: {start!r}-{end!r}")
    return start_idx, end_idx


_CODON = {a + b + c: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[16 * i + 4 * j + k]
          for i, a in enumerate("TCAG") for j, b in enumerate("TCAG") for k, c in enumerate("TCAG")}


def translate(cds: str) -> str:
    return "".join(_CODON.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 2, 3))


def _encodes(cds: str, protein: str) -> bool:
    """CDS encodes exactly `protein` (optionally followed by one stop codon)."""
    L = len(protein)
    if len(cds) not in (3 * L, 3 * L + 3):
        return False
    return translate(cds[:3 * L]) == protein


def _pick_col(columns, preferred: list[str]) -> str | None:
    low = {str(c).strip().lower(): c for c in columns}
    return next((low[p] for p in preferred if p in low), None)


def _read_tf_table(names: list[str], dirs: list[Path], notes: list[str]):
    """Read the first existing table; normalize to columns _id, _cds[, _family]."""
    for folder in dirs:
        for n in names:
            path = folder / n
            if not path.exists():
                continue
            try:
                tbl = pd.read_excel(path) if path.suffix == ".xlsx" else pd.read_csv(path)
            except ImportError:  # .xlsx without openpyxl: try the next candidate (e.g. .csv)
                notes.append(f"{path.name}: install openpyxl to read .xlsx; skipped")
                continue
            id_col, cds_col = _pick_col(tbl.columns, _ID_COLS), _pick_col(tbl.columns, _CDS_COLS)
            if id_col is None or cds_col is None:
                notes.append(f"{path.name}: no ID/CDS columns found in {list(tbl.columns)}")
                continue
            out = pd.DataFrame({
                "_id": tbl[id_col].astype(str).str.strip(),
                "_cds": tbl[cds_col].astype(str).str.upper().str.replace(r"\s+", "", regex=True),
            })
            fam_col = _pick_col(tbl.columns, ["family"])
            if fam_col is not None:
                out["_family"] = tbl[fam_col]
            return out.drop_duplicates("_id"), path
    return None, None


def _empty_pfam() -> pd.DataFrame:
    return pd.DataFrame(columns=["pfam_name", "pfam_acc", "start", "end", "i_evalue", "is_dbd"])


def merge_pfam(hits: pd.DataFrame) -> pd.DataFrame:
    """Merge overlapping or abutting hits that belong to the same domain group.

    Group = PFAM_GROUPS entry (e.g. all MYB models -> "MYB DBD"), otherwise the
    Pfam model name itself (so overlapping copies of one model also collapse).
    Output keeps the hit-table columns; `pfam_name` becomes the group label,
    `pfam_acc`/`members` list the contributing models, `is_dbd` is any().
    """
    cols = ["pfam_name", "pfam_acc", "start", "end", "i_evalue", "is_dbd", "members", "n_hits"]
    if hits.empty:
        return pd.DataFrame(columns=cols)
    h = hits.assign(group=hits.pfam_name.map(lambda n: _PFAM_GROUP_OF.get(n, n)))
    out = []
    for grp, g in h.groupby("group", sort=False):
        g = g.sort_values("start")
        cur = None
        for r in g.itertuples(index=False):
            if cur is not None and r.start <= cur["end"] + 1:
                cur["end"] = max(cur["end"], r.end)
                cur["rows"].append(r)
            else:
                if cur is not None:
                    out.append(cur)
                cur = {"group": grp, "start": r.start, "end": r.end, "rows": [r]}
        out.append(cur)
    rows = []
    for c in out:
        names = list(dict.fromkeys(r.pfam_name for r in c["rows"]))
        accs = list(dict.fromkeys(r.pfam_acc for r in c["rows"]))
        rows.append({
            "pfam_name": c["group"],
            "pfam_acc": ", ".join(accs),
            "start": int(c["start"]),
            "end": int(c["end"]),
            "i_evalue": float(min(r.i_evalue for r in c["rows"])),
            "is_dbd": int(any(r.is_dbd for r in c["rows"])),
            "members": ", ".join(names),
            "n_hits": len(c["rows"]),
        })
    return pd.DataFrame(rows, columns=cols).sort_values("start").reset_index(drop=True)


def _place_tiles(seq: str, positions, tile_seqs, activities) -> pd.DataFrame:
    """Locate each 40-aa empirical tile on the protein.

    `position` is a nominal 1-based start on a 10-aa grid. For the C-terminal
    tile the real tile is the last 40 aa, so we locate tiles by sequence match
    (nearest occurrence to the nominal position). Tiles whose sequence is not
    found (screened from a different isoform) keep their nominal position and
    are flagged placement="nominal".
    """
    rows = []
    for p, s, a in zip(positions, tile_seqs, activities):
        p = int(p)
        placement = "exact"
        if seq[p - 1:p - 1 + len(s)] != s:
            hits, i = [], seq.find(s)
            while i >= 0:
                hits.append(i + 1)
                i = seq.find(s, i + 1)
            if hits:
                p = min(hits, key=lambda h: abs(h - p))
                placement = "shifted"
            else:
                placement = "nominal"
        end = min(p + len(s) - 1, len(seq))
        rows.append((p, end, float(a), s, placement))
    return pd.DataFrame(rows, columns=["start", "end", "activity", "tile_seq", "placement"])


def call_ads(preds: np.ndarray, seq: str,
             threshold: float = AD_THRESHOLD, min_len: int = AD_MIN_LEN) -> pd.DataFrame:
    """Contiguous runs where preds >= threshold, at least `min_len` residues long."""
    above = np.concatenate([[False], np.asarray(preds) >= threshold, [False]])
    edges = np.flatnonzero(np.diff(above.astype(np.int8)))
    rows = []
    for s0, e0 in zip(edges[::2], edges[1::2]):  # 0-based half-open [s0, e0)
        if e0 - s0 >= min_len:
            seg = preds[s0:e0]
            rows.append({
                "AD": f"AD{len(rows) + 1}",
                "Start": int(s0 + 1),
                "End": int(e0),
                "Length": int(e0 - s0),
                "Peak Score": round(float(seg.max()), 3),
                "Mean Score": round(float(seg.mean()), 3),
                "Sequence": seq[s0:e0],
            })
    return pd.DataFrame(rows, columns=["AD", "Start", "End", "Length",
                                       "Peak Score", "Mean Score", "Sequence"])


if __name__ == "__main__":
    import time
    t = time.time()
    data = TFData()
    print(f"loaded in {time.time() - t:.1f}s; {len(data.preds)} TFs")
    for gid in ["Potri.001G018400", "AT1G16060.1", "at1g35240"]:
        rec = data.get(gid)
        print(f"\n{rec.species} {rec.gene_id} ({rec.protein_id}) {rec.family} len={rec.length}")
        print(" tiles:", len(rec.tiles), rec.tiles.placement.value_counts().to_dict())
        print(" pfam:", rec.pfam[["pfam_name", "start", "end"]].values.tolist())
        print(rec.ads.drop(columns="Sequence").to_string(index=False) if len(rec.ads) else " no ADs")
    for q in ["Potri.001G258700", "AT5G12870"]:
        og = data.orthogroup(q)
        print(f"\northogroup({q}): Ptr={og.ptr} Ath={og.ath} ratio={og.ratio} aligned={og.aligned} aln_len={og.aln_len}")
    print(data.orthogroup_summary(data.orthogroup("AT5G12870")).to_string(index=False))
    n_ad = sum(len(call_ads(np.asarray(p, float), s)) > 0 for p, s in zip(data.preds.preds, data.preds.seq))
    print(f"\nTFs with >=1 AD: {n_ad}/{len(data.preds)}")
