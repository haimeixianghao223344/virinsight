"""免疫影响判定：把氨基酸突变映射到表位 / 关键残基 / 结构域 / 糖基化位点。

输出的"免疫影响"为知识库规则推断（evidence=curated 位点 + heuristic 规则），
不等同于实验测定的逃逸。
"""
from __future__ import annotations

# 氨基酸性质（用于启发式严重性评分）
HYDROPHOBIC = set("AVLIMFWP")
POLAR = set("STNQYC")
POS = set("KRH")
NEG = set("DE")
GLY = "G"

BLOSUM62_SIMPLE = {
    # 极简替代打分：命中保守替换得 1，性质跨类得 0（原型；完整版可导入BLOSUM62矩阵）
}


def residue_class(a: str) -> str:
    if a in HYDROPHOBIC:
        return "疏水"
    if a in POS:
        return "正电"
    if a in NEG:
        return "负电"
    if a in POLAR:
        return "极性"
    if a == "G":
        return "甘氨酸(柔性)"
    return "其他"


def severity_score(ref_aa: str, alt_aa: str) -> int:
    """0-3 启发式结构严重性：0保守 / 1中等 / 2较剧烈 / 3终止/缺失/插入。"""
    if alt_aa == "*" or ref_aa == "*":
        return 3
    if alt_aa == "-" or ref_aa == "-":
        return 3
    if ref_aa == alt_aa:
        return 0
    if residue_class(ref_aa) == residue_class(alt_aa):
        return 0
    # 脯氨酸/甘氨酸引入或去除（破坏二级结构）
    if (ref_aa in "PG") != (alt_aa in "PG"):
        return 2
    classes = {residue_class(ref_aa), residue_class(alt_aa)}
    if ("正电" in classes or "负电" in classes) and "疏水" in classes:
        return 2
    if ("正电" in classes and "负电" in classes):
        return 2
    return 1


SEVERITY_LABEL = {0: "保守替换", 1: "中等改变", 2: "剧烈(电荷/极性/构象)", 3: "终止/缺失/插入"}


def _in_ranges(pos: int, ranges: list[tuple[int, int]]) -> list:
    return [(a, b) for a, b in ranges if a <= pos <= b]


def analyze_mutation(protein: str, mut: dict, epi_kb: dict, gene_def: dict) -> dict:
    """对单条氨基酸突变做免疫/结构注释。返回影响标签列表。"""
    tags = []
    evid = []
    pos = mut["aa_pos"]
    alt = mut.get("alt_aa", "")
    mtype = mut["type"]

    pkb = epi_kb.get("proteins", {}).get(protein, {})

    # 1. 关键残基（人工精选，证据强度最高）
    for kr in pkb.get("key_residues", []):
        if isinstance(pos, int) and kr["aa"] == pos:
            tags.append(f"关键功能位点 {kr['label']}：{kr['aspect']}")
            evid.append(("curated", kr.get("effect_note", "")))
            if mtype != "snp" or alt not in ("-", ""):
                tags.append(f"潜在影响：{kr.get('effect_note','')}")

    # 2. 表位区间
    for bucket, label, kind in [
        ("bcell_neutralizing", "中和/B细胞表位", "抗体逃逸"),
        ("tcell_mhc1", "MHC-I T细胞表位", "CD8 T细胞逃逸"),
        ("tcell_mhc2", "MHC-II T细胞表位", "CD4 T细胞应答改变"),
    ]:
        for ep in pkb.get(bucket, []):
            if isinstance(pos, int) and ep["aa_start"] <= pos <= ep["aa_end"]:
                tags.append(f"位于{label}「{ep['name']}」({ep['aa_start']}-{ep['aa_end']})→ 潜在{kind}")
                evid.append(("curated", ep.get("note", "")))

    # 3. 结构域
    for d in gene_def.get("domains", []):
        if isinstance(pos, int) and d["aa_start"] <= pos <= d["aa_end"]:
            tags.append(f"结构域 {d['name']}({d['aa_start']}-{d['aa_end']})：{d['note']}")

    # 4. 糖基化位点 N-X-S/T (X≠P)：获得/丢失
    glyco = set(gene_def.get("glycosylation_sites", []))
    if isinstance(pos, int):
        ref_p = None
        # 已知糖位被改 → 丢失
        if pos in glyco and mut.get("ref_aa") == "N":
            tags.append("丢失已知N-糖基化位点（可能暴露新表位/改变构象）")
        # 检查新糖位获得需要上下文：交给 batch 层面（_glycosylation_gain）

    # 5. furin 切位等特殊位点（通过结构域命名）
    for d in gene_def.get("domains", []):
        if "切位" in d.get("name", "") and isinstance(pos, int) and d["aa_start"] <= pos <= d["aa_end"]:
            tags.append("位于蛋白酶切位附近，可能影响蛋白加工/感染性")

    # 6. 启发式严重性
    if mtype in ("missense", "nonsense") and mut.get("ref_aa") and alt:
        sc = severity_score(mut["ref_aa"], alt if mtype != "nonsense" else "*")
    elif mtype in ("del", "ins"):
        sc = 3
    else:
        sc = 1
    tags.append(f"结构严重性评分 {sc}/3（{SEVERITY_LABEL[sc]}，启发式）")

    return dict(tags=tags, severity=sc, curated_hits=[e for e in evid if e[0] == "curated"])


def glycosylation_gains(ref_prot: str, qry_prot: str) -> list[dict]:
    """检测 N-X-S/T(X≠P) N-糖基化基序的获得与丢失。"""
    import re
    pat = re.compile(r"N[^P][ST]")

    def sites(s):
        return {m.start() + 1 for m in pat.finditer(s)}

    rs, qs = sites(ref_prot), sites(qry_prot)
    gained = sorted(qs - rs)
    lost = sorted(rs - qs)
    return [dict(aa_pos=p, kind="gain", label=f"N{p} 新增N-糖基化基序") for p in gained] + \
           [dict(aa_pos=p, kind="loss", label=f"N{p} 丢失N-糖基化基序") for p in lost]
