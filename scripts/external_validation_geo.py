"""
Prepare GEO cache files used by run_external_validation.py.

This script creates/validates:
  - geo_cache/gse91061_target_expression.csv
  - geo_cache/gse91061_sample_meta.csv

It prefers existing committed processed CSVs for deterministic reproducibility.
If files are missing, it can rebuild them from raw GEO artifacts (and optionally
download missing raw artifacts from NCBI GEO FTP).
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import tarfile
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from typing import Dict

import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(SCRIPT_DIR, "geo_cache")
os.makedirs(CACHE_DIR, exist_ok=True)

TARGET_EXPR_PATH = os.path.join(CACHE_DIR, "gse91061_target_expression.csv")
TARGET_META_PATH = os.path.join(CACHE_DIR, "gse91061_sample_meta.csv")
RAW_FPKM_PATH = os.path.join(CACHE_DIR, "GSE91061_fpkm.csv.gz")
RAW_XML_TGZ_PATH = os.path.join(CACHE_DIR, "GSE91061_family.xml.tgz")
RAW_XML_PATH = os.path.join(CACHE_DIR, "GSE91061_family.xml")

SUPPL_INDEX_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/"
MINIML_XML_TGZ_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/miniml/"
    "GSE91061_family.xml.tgz"
)

# Entrez IDs in GSE91061 FPKM matrix -> symbols expected by run_external_validation.py
ENTREZ_TO_GENE = {
    "29126": "CD274",
    "406": "ARNTL",
    "5187": "PER1",
    "567": "B2M",
    "9575": "CLOCK",
    "8864": "PER2",
    "1407": "CRY1",
    "1408": "CRY2",
    "3105": "HLA_A",
    "3106": "HLA_B",
    "3107": "HLA_C",
    "2697": "GJA1",
    "2706": "GJB2",
    "4609": "MYC",
    "7157": "TP53",
}

REQUIRED_META_COLS = {"sample_id", "patient", "timepoint", "response"}
REQUIRED_EXPR_COLS = set(ENTREZ_TO_GENE.values())


def _download_file(url: str, dst_path: str) -> None:
    print(f"Downloading: {url}")
    with urllib.request.urlopen(url, timeout=120) as resp, open(dst_path, "wb") as out:
        shutil.copyfileobj(resp, out)
    print(f"  Saved: {dst_path}")


def _find_fpkm_url() -> str:
    with urllib.request.urlopen(SUPPL_INDEX_URL, timeout=120) as resp:
        html = resp.read().decode("utf-8", errors="ignore")

    candidates = re.findall(r'href="([^"]*fpkm\.csv\.gz)"', html, flags=re.IGNORECASE)
    if not candidates:
        raise RuntimeError("Could not find FPKM file in GEO supplementary index.")

    # Prefer the known hg19KnownGene fpkm file when present.
    preferred = [c for c in candidates if "KnownGene.fpkm.csv.gz" in c]
    rel_path = preferred[0] if preferred else candidates[0]
    return urllib.parse.urljoin(SUPPL_INDEX_URL, rel_path)


def _normalize_response(value: str) -> str:
    v = (value or "").strip().upper().replace(" ", "")
    if v in {"PRCR", "PR/CR", "PR", "CR"}:
        return "PRCR"
    if v == "PD":
        return "PD"
    if v == "SD":
        return "SD"
    return value.strip() if value else ""


def _extract_xml_from_tgz(tgz_path: str, xml_path: str) -> None:
    with tarfile.open(tgz_path, "r:gz") as tar:
        xml_members = [m for m in tar.getmembers() if m.name.lower().endswith(".xml")]
        if not xml_members:
            raise RuntimeError(f"No XML file found inside {tgz_path}.")
        member = xml_members[0]
        with tar.extractfile(member) as src, open(xml_path, "wb") as dst:
            if src is None:
                raise RuntimeError(f"Failed to extract {member.name} from {tgz_path}.")
            shutil.copyfileobj(src, dst)
    print(f"  Extracted XML: {xml_path}")


def _load_response_map(xml_path: str) -> Dict[str, str]:
    ns = {"geo": "http://www.ncbi.nlm.nih.gov/geo/info/MINiML"}
    tree = ET.parse(xml_path)
    root = tree.getroot()

    response_map: Dict[str, str] = {}
    for sample in root.findall(".//geo:Sample", ns):
        title = sample.findtext("geo:Title", "", ns).strip()
        response = ""
        for ch in sample.findall(".//geo:Characteristics", ns):
            if (ch.get("tag") or "").strip().lower() == "response":
                response = _normalize_response(ch.text or "")
                break
        if title:
            response_map[title] = response
    return response_map


def _build_target_expression(raw_fpkm_path: str) -> pd.DataFrame:
    raw = pd.read_csv(raw_fpkm_path, index_col=0)
    raw.index = raw.index.astype(str)

    missing = [eid for eid in ENTREZ_TO_GENE if eid not in raw.index]
    if missing:
        raise RuntimeError(f"Missing expected Entrez IDs in FPKM matrix: {missing}")

    expr = raw.loc[list(ENTREZ_TO_GENE.keys())].copy()
    expr.index = [ENTREZ_TO_GENE[eid] for eid in expr.index]
    expr = expr.T
    expr.index.name = "sample_id"
    return expr


def _build_sample_meta(sample_ids: pd.Index, response_map: Dict[str, str]) -> pd.DataFrame:
    records = []
    for sid in sample_ids:
        parts = str(sid).split("_")
        patient = parts[0] if len(parts) > 0 else ""
        timepoint = parts[1] if len(parts) > 1 else ""
        records.append(
            {
                "sample_id": sid,
                "patient": patient,
                "timepoint": timepoint,
                "response": _normalize_response(response_map.get(str(sid), "")),
            }
        )
    meta = pd.DataFrame(records).set_index("sample_id")
    return meta


def _validate_processed_outputs() -> bool:
    if not (os.path.exists(TARGET_EXPR_PATH) and os.path.exists(TARGET_META_PATH)):
        return False

    expr = pd.read_csv(TARGET_EXPR_PATH, index_col=0)
    meta = pd.read_csv(TARGET_META_PATH, index_col=0)

    expr_ok = REQUIRED_EXPR_COLS.issubset(set(expr.columns))
    meta_ok = REQUIRED_META_COLS.issubset(set(["sample_id"] + list(meta.columns)))
    overlap_ok = len(set(expr.index).intersection(set(meta.index))) > 0
    return expr_ok and meta_ok and overlap_ok


def prepare_gse91061_cache(force_refresh: bool = False, allow_download: bool = True) -> None:
    """Ensure processed GSE91061 cache files exist and are consistent."""
    if not force_refresh and _validate_processed_outputs():
        print("Processed GEO cache already present and valid.")
        print(f"  {TARGET_EXPR_PATH}")
        print(f"  {TARGET_META_PATH}")
        return

    if force_refresh:
        print("Force refresh requested; rebuilding processed GEO cache files.")
    else:
        print("Processed GEO cache missing/incomplete; attempting rebuild.")

    if force_refresh or not os.path.exists(RAW_FPKM_PATH):
        if not allow_download:
            raise RuntimeError(
                f"Raw FPKM file missing and downloads disabled: {RAW_FPKM_PATH}"
            )
        fpkm_url = _find_fpkm_url()
        _download_file(fpkm_url, RAW_FPKM_PATH)

    if force_refresh or not os.path.exists(RAW_XML_PATH):
        if force_refresh or not os.path.exists(RAW_XML_TGZ_PATH):
            if not allow_download:
                raise RuntimeError(
                    f"Raw MINiML tgz missing and downloads disabled: {RAW_XML_TGZ_PATH}"
                )
            _download_file(MINIML_XML_TGZ_URL, RAW_XML_TGZ_PATH)
        _extract_xml_from_tgz(RAW_XML_TGZ_PATH, RAW_XML_PATH)

    expr = _build_target_expression(RAW_FPKM_PATH)
    response_map = _load_response_map(RAW_XML_PATH)
    meta = _build_sample_meta(expr.index, response_map)

    expr.to_csv(TARGET_EXPR_PATH)
    meta.to_csv(TARGET_META_PATH)

    print(f"Saved: {TARGET_EXPR_PATH} (rows={len(expr)}, cols={len(expr.columns)})")
    print(f"Saved: {TARGET_META_PATH} (rows={len(meta)}, cols={len(meta.columns)})")

    if not _validate_processed_outputs():
        raise RuntimeError("Processed GEO cache files were written but failed validation.")

    print("GSE91061 cache preparation complete.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare GSE91061 cache files.")
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Rebuild processed files from raw assets (re-download raw files when needed).",
    )
    parser.add_argument(
        "--no-download",
        action="store_true",
        help="Do not download raw GEO files; only use local cache files.",
    )
    args = parser.parse_args()

    print("=" * 65)
    print("  GEO Cache Prep: GSE91061 (Riaz et al. 2017)")
    print("=" * 65)

    try:
        prepare_gse91061_cache(
            force_refresh=args.refresh,
            allow_download=not args.no_download,
        )
    except Exception as exc:
        print(f"Cache preparation failed: {exc}")
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
