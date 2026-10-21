"""数据与知识库加载。"""
from __future__ import annotations
import json
from pathlib import Path
from Bio import SeqIO

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
REF_DIR = DATA_DIR / "refs"
KB_DIR = DATA_DIR / "kb"


def load_fasta(path: str | Path):
    """返回 [(record_id, description, sequence_upper), ...]。"""
    out = []
    for r in SeqIO.parse(str(path), "fasta"):
        out.append((r.id, r.description, str(r.seq).upper()))
    if not out:
        raise ValueError(f"FASTA 中未读到序列: {path}")
    return out


def parse_fasta_text(text: str):
    """从粘贴的 FASTA 文本解析序列。"""
    from io import StringIO
    out = []
    for r in SeqIO.parse(StringIO(text.strip()), "fasta"):
        out.append((r.id, r.description, str(r.seq).upper()))
    return out


def load_json(name: str):
    with open(KB_DIR / name, encoding="utf-8") as f:
        return json.load(f)


def load_viruses():
    return load_json("viruses.json")


def get_virus(vid: str):
    for v in load_viruses():
        if v["id"] == vid:
            return v
    raise KeyError(f"未知病毒: {vid}")


def load_reference(vid: str) -> str:
    v = get_virus(vid)
    recs = load_fasta(REF_DIR / v["ref_file"])
    return recs[0][2]
