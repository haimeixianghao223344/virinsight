"""端到端分析流水线：识别 → 变异 → 注释 → 免疫 → 结构影响 → 谱系签名。"""
from __future__ import annotations
from .io import load_reference, get_virus, load_json
from .identify import identify
from .variant import call_variants, summarize
from .annotate import annotate_genes, count_syn_nonsyn
from .immuno import analyze_mutation, glycosylation_gains
from .online import annotate_online


def identify_virus(seq: str):
    cands = identify(seq)
    return cands


def analyze_sample(ref_seq: str, qry_seq: str, virus_def: dict,
                   fetch_online: bool = False, online_cache: dict | None = None) -> dict:
    genes = load_json(virus_def["genes_kb"])
    epi = load_json(virus_def["epitopes_kb"])
    variants = call_variants(ref_seq, qry_seq)
    gstats = summarize(variants)
    annotated = annotate_genes(ref_seq, qry_seq, genes, variants)
    snds = count_syn_nonsyn(ref_seq, qry_seq, genes)

    # 组装每基因的蛋白突变 + 免疫/结构影响
    gene_rows = []
    impact_rows = []
    for g in genes:
        a = annotated[g["name"]]
        muts = a["mutations"]
        # 糖基化获得/丢失
        gly = glycosylation_gains(a["ref_protein"], a["qry_protein"])
        rows = []
        for m in muts:
            im = analyze_mutation(g["name"], m, epi, g)
            rows.append(dict(mutation=m, impact=im))
            impact_rows.append(dict(
                gene=g["name"], protein=g["name"], label=m["label"],
                aa_pos=m["aa_pos"], ref_aa=m.get("ref_aa"), alt_aa=m.get("alt_aa"),
                type=m["type"], severity=im["severity"],
                immune_tags="; ".join(t for t in im["tags"]),
                in_epitope=any("表位" in t for t in im["tags"]),
                in_key_site=bool(im["curated_hits"]),
            ))
        # 糖位事件也加入影响行
        for gm in gly:
            impact_rows.append(dict(
                gene=g["name"], protein=g["name"], label=gm["label"],
                aa_pos=gm["aa_pos"], ref_aa="N", alt_aa="?",
                type="glyco_" + gm["kind"], severity=2,
                immune_tags=gm["label"] + "（糖基化改变可遮蔽/暴露表位并影响折叠，启发式）",
                in_epitope=False, in_key_site=False,
            ))
        n_epi = sum(1 for r in rows if r["impact"]["curated_hits"] or
                    any("表位" in t for t in r["impact"]["tags"]))
        gene_rows.append(dict(
            name=g["name"], category=g["category"],
            immune_relevance=g["immune_relevance"],
            function=g["function"], immuno_note=g.get("immuno_note", ""),
            n_missense=a["n_missense"], n_indel=a["n_indel"],
            n_immuno_hits=n_epi,
            snds=snds.get(g["name"], {}),
            mutations=[r["mutation"]["label"] for r in rows],
        ))

    # 免疫主作用基因排序：免疫相关度 + 是否实际命中表位/关键位点
    ranked = sorted(gene_rows,
                    key=lambda r: (r["n_immuno_hits"] * 10 + r["immune_relevance"],
                                   r["n_missense"]), reverse=True)

    # 谱系签名匹配（如有）
    lineage_matches = []
    for lin in epi.get("lineage_signatures", []):
        hit = 0; total = 0; detail = []
        for spec in lin["mutations"]:
            total += 1
            ann = annotated.get(spec["protein"])
            if not ann:
                continue
            want_pos, want_alt = spec["aa"], spec["alt"]
            found = None
            for m in ann["mutations"]:
                if m["aa_pos"] == want_pos:
                    found = m
                    break
            ok = False
            if found:
                if want_alt == "del":
                    ok = found["type"] == "del"
                elif want_alt is None:
                    ok = True  # 该位点存在任意突变即计入弱特征
                else:
                    ok = found.get("alt_aa") == want_alt
            if ok:
                hit += 1
                detail.append(f"{spec['protein']}:{want_pos}{want_alt or 'X'}")
        if hit >= 3 or (total and hit / max(total, 1) >= 0.6):
            lineage_matches.append(dict(lineage=lin["lineage"], hit=hit, total=total,
                                        fraction=round(hit / max(total, 1), 2),
                                        matched=detail, note=lin.get("note", "")))

    # 在线增强
    online = {}
    if fetch_online:
        for prot, acc in virus_def.get("uniprot_focus", {}).items():
            ann = annotated.get(prot)
            if not ann:
                continue
            # 只为有突变的重点蛋白查询，控制请求
            if ann["n_missense"] + ann["n_indel"] > 0 or prot in ("S", "Env_gp160"):
                if online_cache is not None and acc in online_cache:
                    online[prot] = online_cache[acc]
                else:
                    online[prot] = annotate_online(acc)
                    if online_cache is not None:
                        online_cache[acc] = online[prot]

    return dict(
        n_variants=len(variants), variant_summary=gstats,
        variants=variants,
        gene_rows=gene_rows, ranked_genes=ranked,
        impact_rows=sorted(impact_rows, key=lambda x: (not x["in_key_site"],
                                                        not x["in_epitope"],
                                                        -x["severity"])),
        lineage_matches=sorted(lineage_matches, key=lambda x: -x["fraction"]),
        online=online,
        annotated=annotated,
    )


def analyze_many(ref_seq: str, samples: list[tuple[str, str]], virus_def: dict) -> dict:
    """多样本(不同时期)分析。samples: [(name, seq), ...]，按传入顺序视为时间序。"""
    per_sample = []
    for name, seq in samples:
        r = analyze_sample(ref_seq, seq, virus_def)
        per_sample.append((name, r))

    # 突变时间线：基因:蛋白突变标签 -> 各样本是否存在
    timeline = {}
    for name, r in per_sample:
        for row in r["impact_rows"]:
            key = (row["gene"], row["aa_pos"], row["ref_aa"], row["alt_aa"], row["type"])
            timeline.setdefault(key, {"gene": row["gene"], "aa_pos": row["aa_pos"],
                                      "ref_aa": row["ref_aa"], "alt_aa": row["alt_aa"],
                                      "type": row["type"], "label": row["label"],
                                      "samples": {}})
            timeline[key]["samples"][name] = 1

    names = [n for n, _ in per_sample]
    timeline_rows = []
    for _, t in timeline.items():
        row = dict(gene=t["gene"], label=t["label"], aa_pos=t["aa_pos"],
                   present=[t["samples"].get(n, 0) for n in names])
        first = next((i for i, p in enumerate(row["present"]) if p), None)
        last = len(row["present"]) - 1 - next((i for i, p in enumerate(reversed(row["present"])) if p), 0)
        row["first_seen"] = names[first] if first is not None else ""
        row["last_seen"] = names[last] if last is not None else ""
        row["fixed"] = row["present"][-1] == 1 if row["present"] else 0
        timeline_rows.append(row)
    timeline_rows.sort(key=lambda r: (r["gene"], r["aa_pos"]))

    return dict(sample_names=names, per_sample=per_sample, timeline_rows=timeline_rows)
