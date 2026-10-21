"""基因注释：把核酸变异映射到氨基酸突变，并估算非同义/同义突变。"""
from __future__ import annotations
from Bio.Seq import Seq
from difflib import SequenceMatcher
from .variant import call_variants, reverse_complement

CODON_TABLE = 1  # 标准遗传密码（原型；冠状病毒移码/特殊密码子在README中说明局限）

AA3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
       "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
       "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
       "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
       "*": "Ter", "X": "X", "-": "-"}


def _apply_variants_to_segment(seg_ref: str, seg_start0: int, variants: list[dict]) -> str:
    """按参考坐标把变异施加到一段参考子序列上（0-based 全局起点 seg_start0）。"""
    edits = []  # (局部起点, 局部终点, 替换文本)
    for v in variants:
        a, b = v["ref_start"], v["ref_end"]
        if b <= seg_start0 or a >= seg_start0 + len(seg_ref):
            continue
        la, lb = max(a, seg_start0) - seg_start0, min(b, seg_start0 + len(seg_ref)) - seg_start0
        if v["type"] == "snp":
            if la == lb:  # 插入落在边界，跳过边界外
                if not (a == b and seg_start0 <= a <= seg_start0 + len(seg_ref)):
                    continue
            edits.append((la, lb, v["alt_nt"]))
        elif v["type"] == "del":
            edits.append((la, lb, ""))
        elif v["type"] == "ins":
            edits.append((la, lb, v["alt_nt"]))
    edits.sort(key=lambda e: (e[0], -e[1]))
    out, cur = [], 0
    for a, b, rep in edits:
        if a < cur:  # 重叠编辑（罕见复杂indel），跳过后者
            continue
        out.append(seg_ref[cur:a]); out.append(rep); cur = b
    out.append(seg_ref[cur:])
    return "".join(out)


def translate(seq: str, phase: int = 1) -> str:
    """phase: 基因起始在密码子中的相位(1/2/3)。原型按完整读码框从基因起点翻译。"""
    s = seq[(phase - 1):]
    return str(Seq(s).translate(table=CODON_TABLE, to_stop=False))


def protein_mutations(ref_prot: str, qry_prot: str) -> list[dict]:
    """比对两蛋白，输出氨基酸突变列表（位置相对参考蛋白，1-based）。"""
    sm = SequenceMatcher(None, ref_prot, qry_prot, autojunk=False)
    out = []
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == "equal":
            continue
        rp, qp = ref_prot[i1:i2], qry_prot[j1:j2]
        if tag == "replace" and len(rp) == len(qp):
            for off, (a, b) in enumerate(zip(rp, qp)):
                pos = i1 + off + 1
                out.append(dict(aa_pos=pos, ref_aa=a, alt_aa=b,
                                label=f"{AA3.get(a,a)}{pos}{AA3.get(b,b)}",
                                type="missense" if b != "*" else "nonsense"))
        elif tag == "delete":
            out.append(dict(aa_pos=i1 + 1, ref_aa=rp, alt_aa="-",
                            label=f"Δ{len(rp)}aa@{i1+1}", type="del",
                            detail=f"{AA3.get(rp[0],rp[0])}{i1+1}-{AA3.get(rp[-1],rp[-1])}{i2}缺失"))
        elif tag == "insert":
            out.append(dict(aa_pos=i1, ref_aa="-", alt_aa=qp,
                            label=f"ins{len(qp)}aa@{i1}", type="ins"))
        else:
            # 不等长替换：按短端逐点 + 额外indel
            n = min(len(rp), len(qp))
            for off in range(n):
                if rp[off] != qp[off]:
                    pos = i1 + off + 1
                    out.append(dict(aa_pos=pos, ref_aa=rp[off], alt_aa=qp[off],
                                    label=f"{AA3.get(rp[off],rp[off])}{pos}{AA3.get(qp[off],qp[off])}",
                                    type="missense"))
            if len(rp) > n:
                out.append(dict(aa_pos=i1 + n + 1, ref_aa=rp[n:], alt_aa="-",
                                label=f"Δ{len(rp)-n}aa@{i1+n+1}", type="del"))
            else:
                out.append(dict(aa_pos=i1 + n, ref_aa="-", alt_aa=qp[n:],
                                label=f"ins{len(qp)-n}aa@{i1+n}", type="ins"))
    out.sort(key=lambda x: x["aa_pos"])
    return out


def annotate_genes(ref: str, qry: str, genes: list[dict], variants: list[dict]) -> dict:
    """对每个基因提取参考/查询片段，翻译，调用氨基酸突变。返回 gene -> 结果 dict。"""
    result = {}
    for g in genes:
        s0, s1 = g["nt_start"] - 1, g["nt_end"]
        strand = g.get("strand", 1)
        phase = g.get("phase", 1)
        ref_seg = ref[s0:s1] if strand == 1 else reverse_complement(ref[s0:s1])
        qry_seg = qry[s0:s1] if strand == 1 else reverse_complement(qry[s0:s1])
        if strand == 1:
            qry_seg = _apply_variants_to_segment(ref_seg, s0, variants)
        # 注：负链基因的变异施加为原型简化（内置库均为正链），直接用qry提取片段
        ref_p = translate(ref_seg, phase)
        qry_p = translate(qry_seg, phase)
        pm = protein_mutations(ref_p.rstrip("*"), qry_p.rstrip("*"))
        result[g["name"]] = dict(
            gene=g, ref_protein=ref_p, qry_protein=qry_p,
            mutations=pm,
            n_missense=sum(1 for m in pm if m["type"] in ("missense", "nonsense")),
            n_indel=sum(1 for m in pm if m["type"] in ("del", "ins")),
        )
    return result


def count_syn_nonsyn(ref: str, qry: str, genes: list[dict]) -> dict:
    """用 dNdScv 思想的极简化版：仅统计基因内 SNP 导致的同义和非同义密码子改变。
    返回 {gene: {syn, nonsyn, stop}}；用于 dN/dS 粗估（未做位点标准化）。"""
    out = {}
    variants = call_variants(ref, qry)
    for g in genes:
        if g.get("strand", 1) != 1:
            continue
        s0, s1 = g["nt_start"] - 1, g["nt_end"]
        ph = g.get("phase", 1) - 1
        syn = nonsyn = stop = 0
        # 按密码子收集SNP
        codon_snps: dict[int, list[tuple[int, str, str]]] = {}
        for v in variants:
            if v["type"] != "snp":
                continue
            pos0 = v["ref_pos"] - 1
            if not (s0 <= pos0 < s1):
                continue
            frame = (pos0 - s0 - ph) % 3
            if (pos0 - s0 - ph) < 0:
                continue
            cidx = (pos0 - s0 - ph) // 3
            codon_snps.setdefault(cidx, []).append((frame, v["ref_nt"], v["alt_nt"]))
        for cidx, snps in codon_snps.items():
            cstart = s0 + ph + cidx * 3
            if cstart + 3 > s1:
                continue
            rc = ref[cstart:cstart + 3]
            qc = list(rc)
            for frame, a, b in snps:
                if qc[frame] == a:
                    qc[frame] = b
            qc = "".join(qc)
            if "N" in qc:
                continue
            ra = str(Seq(rc).translate())
            qa = str(Seq(qc).translate())
            if ra == qa:
                syn += 1
            elif qa == "*":
                stop += 1
            else:
                nonsyn += 1
        out[g["name"]] = dict(syn=syn, nonsyn=nonsyn, stop=stop)
    return out
