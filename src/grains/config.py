"""Paths and settings shared across the pipeline."""

from __future__ import annotations

from pathlib import Path

# --- paths -----------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "data"
REFERENCE = DATA / "reference"
RAW = DATA / "raw"
PROCESSED = DATA / "processed"
RESULTS = REPO_ROOT / "results"

GRAINS_CSV = REFERENCE / "grains.csv"
PEPTIDES_CSV = REFERENCE / "peptides.csv"

# --- UniProt -------------------------------------------------------------------
UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_PAGE_SIZE = 200          # max 500; kept modest to be gentle on the API
REQUEST_TIMEOUT = 30            # seconds
CONTACT_EMAIL = "contact.s.foster@gmail.com"  # sent as User-Agent per UniProt etiquette

# Seed-storage protein search terms per family. Keys match `storage_protein_family`
# tokens in grains.csv (first word, lowercased). Values are OR-ed protein-name queries.
STORAGE_PROTEIN_TERMS: dict[str, list[str]] = {
    "gliadin/glutenin": ["gliadin", "glutenin", "avenin-like", "gluten"],
    "hordein": ["hordein"],
    "secalin": ["secalin"],
    "avenin": ["avenin"],
    "glutelin/orzenin": ["glutelin", "prolamin", "orzenin", "seed storage"],
    "zein": ["zein"],
    "kafirin": ["kafirin"],
    "setarin": ["prolamin", "seed storage", "setarin"],
    "panicin": ["prolamin", "seed storage", "panicin"],
    "prolamin": ["prolamin", "seed storage"],
    "13s": ["globulin", "2S albumin", "13S", "vicilin", "legumin", "seed storage"],
    "11s/2s": ["globulin", "2S albumin", "11S", "legumin", "seed storage"],
    "11s": ["globulin", "11S", "legumin", "seed storage"],
}

# Only pull entries whose sequence length falls in this window (drop obvious
# fragments and unrelated large hits). Seed prolamins are typically ~100-500 aa.
MIN_SEQ_LEN = 60
MAX_SEQ_LEN = 900


def ensure_dirs() -> None:
    for d in (RAW, PROCESSED, RESULTS):
        d.mkdir(parents=True, exist_ok=True)
