"""病毒种类识别：基于参考库 k-mer 轮廓。"""
from __future__ import annotations
from .io import load_viruses, load_fasta, REF_DIR


def _kmers(seq: str, k: int = 21) -> set[str]:
    seq = seq.upper()
    if len(seq) < k:
        return {seq}
    return {seq[i:i + k] for i in range(len(seq) - k + 1) if "N" not in seq[i:i + k]}


def identify(seq: str, k: int = 21, min_identity: float = 0.15):
    """返回 [(virus_id, name, jaccard, shared_fraction_of_sample), ...]，按相似度降序。"""
    sk = _kmers(seq, k)
    results = []
    for v in load_viruses():
        ref = load_fasta(REF_DIR / v["ref_file"])[0][2]
        rk = _kmers(ref, k)
        if not sk or not rk:
            continue
        inter = len(sk & rk)
        jac = inter / len(sk | rk)
        frac = inter / len(sk)
        results.append((v["id"], v["name"], jac, frac, v))
    results.sort(key=lambda x: (x[3], x[2]), reverse=True)
    return results


def best_match(seq: str):
    r = identify(seq)
    if not r:
        return None
    return r[0]
