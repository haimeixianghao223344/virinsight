"""近缘基因组快速变异调用（核酸层）。

使用 edlib（C 扩展，edit-distance + 回溯 CIGAR）对 30kb 病毒基因组做全局比对，
毫秒级得到带 gap 的精确比对，再逐列解析为 SNP/缺失/插入。
定位：病毒株内/株系间追踪；不处理大段重排。
"""
from __future__ import annotations
import re
import edlib
from Bio.Align import PairwiseAligner

_COMP = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _gapped(ref: str, qry: str):
    res = edlib.align(qry, ref, mode="NW", task="path")
    # edlib: query 为第一参数，target 为第二；cigar 沿 query
    cigar = res.get("cigar")
    if cigar is None:
        raise RuntimeError("edlib 未返回比对路径")
    parsed = [(int(n), op) for n, op in re.findall(r"(\d+)([=XDI])", cigar)]
    ri = qi = 0
    ref_cols: list[int | None] = []
    qry_cols: list[int | None] = []
    for length, op in parsed:
        for _ in range(length):
            if op == "=":      # match
                ref_cols.append(ri); qry_cols.append(qi); ri += 1; qi += 1
            elif op == "X":    # mismatch
                ref_cols.append(ri); qry_cols.append(qi); ri += 1; qi += 1
            elif op == "D":    # deletion relative to query => ref has base
                ref_cols.append(ri); qry_cols.append(None); ri += 1
            elif op == "I":    # insertion relative to query => query has base
                ref_cols.append(None); qry_cols.append(qi); qi += 1
    return ref_cols, qry_cols


_refine = PairwiseAligner()
_refine.mode = "global"
_refine.match_score = 2
_refine.mismatch_score = -2
_refine.open_gap_score = -6
_refine.extend_gap_score = -3


def _block_events(rseg: str, qseg: str, ref_offset: int) -> list[dict]:
    """对单个差异小块用高 gap 罚分重比对（偏好错配而非相邻 indel），展开为事件。"""
    aln = _refine.align(rseg, qseg)[0]
    ti, qi = aln.indices
    out = []
    n = len(ti)
    k = 0
    while k < n:
        rc, qc = ti[k], qi[k]
        if rc != -1 and qc != -1 and rseg[rc] != qseg[qc]:
            j = k
            while j + 1 < n and ti[j + 1] != -1 and qi[j + 1] != -1 \
                    and rseg[ti[j + 1]] != qseg[qi[j + 1]]:
                j += 1
            for t in range(k, j + 1):
                rp = ref_offset + ti[t]
                out.append(dict(type="snp", ref_pos=rp + 1, ref_nt=rseg[ti[t]],
                                alt_nt=qseg[qi[t]], ref_start=rp, ref_end=rp + 1))
            k = j + 1
        elif rc != -1 and qc == -1:
            j = k
            while j + 1 < n and qi[j + 1] == -1:
                j += 1
            rs = ref_offset + rc
            out.append(dict(type="del", ref_pos=rs + 1,
                            ref_nt="".join(rseg[ti[t]] for t in range(k, j + 1)),
                            alt_nt="-", ref_start=rs, ref_end=rs + (j - k) + 1))
            k = j + 1
        elif rc == -1 and qc != -1:
            j = k
            while j + 1 < n and ti[j + 1] == -1:
                j += 1
            prev = (ref_offset + ti[k - 1] + 1) if k > 0 and ti[k - 1] != -1 else 0
            out.append(dict(type="ins", ref_pos=prev, ref_nt="-",
                            alt_nt="".join(qseg[qi[t]] for t in range(k, j + 1)),
                            ref_start=prev, ref_end=prev))
            k = j + 1
        else:
            k += 1
    return out


def call_variants(ref: str, qry: str) -> list[dict]:
    ref = ref.upper().replace("U", "T")
    qry = qry.upper().replace("U", "T")
    ref_cols, qry_cols = _gapped(ref, qry)
    n = len(ref_cols)

    def is_bad(t):
        rc, qc = ref_cols[t], qry_cols[t]
        return rc is None or qc is None or ref[rc] != qry[qc]

    bad = [is_bad(t) for t in range(n)]
    variants = []
    k = 0
    while k < n:
        if not bad[k]:
            k += 1
            continue
        # 合并差异区间，吞并中间 <=2 列的 equal
        j = k
        while j + 1 < n:
            if bad[j + 1]:
                j += 1
                continue
            gap_eq = 0
            t = j + 1
            while t < n and not bad[t]:
                gap_eq += 1
                t += 1
            if t < n and gap_eq <= 2:
                j = t
                continue
            break
        rpos = [ref_cols[t] for t in range(k, j + 1) if ref_cols[t] is not None]
        qpos = [qry_cols[t] for t in range(k, j + 1) if qry_cols[t] is not None]
        if rpos and qpos:
            rseg = ref[rpos[0]:rpos[-1] + 1]
            qseg = qry[qpos[0]:qpos[-1] + 1]
            variants.extend(_block_events(rseg, qseg, rpos[0]))
        elif rpos:
            variants.append(dict(type="del", ref_pos=rpos[0] + 1,
                                 ref_nt="".join(ref[p] for p in rpos), alt_nt="-",
                                 ref_start=rpos[0], ref_end=rpos[-1] + 1))
        else:
            prev = (ref_cols[k - 1] + 1) if k > 0 and ref_cols[k - 1] is not None else 0
            variants.append(dict(type="ins", ref_pos=prev, ref_nt="-",
                                 alt_nt="".join(qry[p] for p in qpos),
                                 ref_start=prev, ref_end=prev))
        k = j + 1
    variants.sort(key=lambda x: (x["ref_pos"], 0 if x["type"] == "snp" else 1))
    return variants


def summarize(variants: list[dict]) -> dict:
    s = {"snp": 0, "del": 0, "ins": 0, "del_nt": 0, "ins_nt": 0}
    for v in variants:
        t = v["type"]
        s[t] = s.get(t, 0) + 1
        if t == "del":
            s["del_nt"] += len(v["ref_nt"])
        elif t == "ins":
            s["ins_nt"] += len(v["alt_nt"])
    return s
