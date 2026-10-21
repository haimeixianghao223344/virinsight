"""HTML 报告导出（自包含，表格+内联样式，无外部依赖）。"""
from __future__ import annotations
import html
import datetime


def _table(headers, rows, cls=""):
    th = "".join(f"<th>{html.escape(str(h))}</th>" for h in headers)
    trs = []
    for r in rows:
        tds = "".join(
            f"<td>{'' if c is None else html.escape(str(c))}</td>" for c in r)
        trs.append(f"<tr>{tds}</tr>")
    return f'<table class="{cls}"><thead><tr>{th}</tr></thead><tbody>{"".join(trs)}</tbody></table>'


def build_single_html(meta: dict, result: dict, virus_def: dict) -> str:
    css = """
    body{font-family:-apple-system,'Segoe UI',Arial,sans-serif;margin:28px;color:#1d2733}
    h1{font-size:22px} h2{font-size:17px;margin-top:26px;border-left:4px solid #2563eb;padding-left:8px}
    table{border-collapse:collapse;width:100%;font-size:12.5px;margin:8px 0}
    th,td{border:1px solid #d4dbe3;padding:5px 8px;text-align:left;vertical-align:top}
    th{background:#eef3fb} .sev3{background:#fde2e1}.sev2{background:#fdf1d6}.key{font-weight:700;color:#b91c1c}
    .muted{color:#64748b;font-size:12px} .tag{display:inline-block;background:#e0ecff;border-radius:4px;padding:1px 6px;margin:1px}
    """
    vs = result["variant_summary"]
    parts = [f"<html><head><meta charset='utf-8'><style>{css}</style></head><body>"]
    parts.append(f"<h1>VirInsight 单基因组免疫分析报告</h1>")
    parts.append(f"<p class='muted'>生成时间 {datetime.datetime.utcnow().isoformat(timespec='seconds')}Z · "
                 f"参考株：{html.escape(virus_def['ref_name'])} · 样本：{html.escape(meta.get('sample',''))} · "
                 f"长度 {meta.get('length','')} nt</p>")
    parts.append(f"<p>核酸变异：SNP <b>{vs['snp']}</b>，缺失 {vs['del']}（{vs['del_nt']} nt），"
                 f"插入 {vs['ins']}（{vs['ins_nt']} nt）</p>")

    if result["lineage_matches"]:
        parts.append("<h2>谱系/关注株签名匹配</h2>")
        parts.append(_table(["谱系", "命中/总位点", "比例", "说明"],
                           [[m["lineage"], f"{m['hit']}/{m['total']}", m["fraction"], m["note"]]
                            for m in result["lineage_matches"]]))

    parts.append("<h2>免疫主作用基因排序（免疫相关度×实际命中表位/功能位点）</h2>")
    rows = []
    for g in result["ranked_genes"]:
        if g["n_missense"] + g["n_indel"] + g["n_immuno_hits"] == 0 and g["immune_relevance"] < 4:
            continue
        rows.append([g["name"], {"structural_antigen": "结构/抗原", "immune_antagonist": "免疫拮抗",
                                 "replication": "复制酶", "structural": "结构", "accessory": "辅助",
                                 "regulatory": "调控"}.get(g["category"], g["category"]),
                     "★" * g["immune_relevance"], g["n_missense"], g["n_indel"], g["n_immuno_hits"],
                     g["immuno_note"]])
    parts.append(_table(["基因", "类别", "免疫相关度", "错义", "蛋白indel", "命中免疫位点", "免疫作用说明"], rows))

    parts.append("<h2>蛋白突变 → 结构/免疫影响明细</h2>")
    rows = []
    for r in result["impact_rows"]:
        cls = "key" if r["in_key_site"] else ""
        rows.append([r["gene"], r["label"], {"missense": "错义", "nonsense": "无义",
                                             "del": "缺失", "ins": "插入"}.get(r["type"], r["type"]),
                     r["severity"], r["immune_tags"]])
    parts.append(_table(["基因/蛋白", "突变", "类型", "结构严重性", "影响判定（证据:知识库精选+启发式）"], rows,
                        "impact"))

    if result.get("online"):
        parts.append("<h2>在线数据库增强注释</h2>")
        for prot, a in result["online"].items():
            u = a["uniprot"]
            if u.get("ok"):
                feats = [[f["type"], f"{f['start']}-{f['end']}", f["description"]] for f in u["features"][:40]]
                parts.append(f"<p><b>{prot}</b> — UniProt {u['accession']} "
                             f"({html.escape(u.get('protein_name',''))})，RCSB 结构: "
                             f"{', '.join(a['rcsb'].get('structures', [])[:10]) or '未检索到'}</p>")
                parts.append(_table(["feature", "位置", "描述"], feats))
            else:
                parts.append(f"<p class='muted'>{prot}: 在线注释不可用（{u.get('error','network')}），已降级为本地分析</p>")
    parts.append("<p class='muted'>免责声明：本报告为科研辅助输出。表位逃逸与结构严重性为知识库规则与在线数据库"
                 "推断，并非ΔΔG计算或实验验证，不构成临床或诊疗结论。通用病毒的分析深度取决于该病毒知识库覆盖。</p>")
    parts.append("</body></html>")
    return "".join(parts)


def build_timeseries_html(meta: dict, many: dict, virus_def: dict) -> str:
    css = """
    body{font-family:-apple-system,'Segoe UI',Arial,sans-serif;margin:28px;color:#1d2733}
    h1{font-size:22px} h2{font-size:17px;margin-top:26px;border-left:4px solid #059669;padding-left:8px}
    table{border-collapse:collapse;font-size:12.5px;margin:8px 0}
    th,td{border:1px solid #d4dbe3;padding:5px 8px;text-align:center}
    th{background:#eafaf3} td.lbl{text-align:left;font-weight:600}
    .on{background:#34a853;color:#fff;border-radius:3px} .off{color:#cbd5e1}
    .muted{color:#64748b;font-size:12px}
    """
    names = many["sample_names"]
    parts = [f"<html><head><meta charset='utf-8'><style>{css}</style></head><body>"]
    parts.append("<h1>VirInsight 多时期基因组突变追踪报告</h1>")
    parts.append(f"<p class='muted'>生成时间 {datetime.datetime.utcnow().isoformat(timespec='seconds')}Z · "
                 f"参考株 {html.escape(virus_def['ref_name'])} · {len(names)} 个时间点："
                 f"{html.escape(' → '.join(names))}</p>")
    head = ["基因", "突变"] + names + ["首次出现", "末期固定"]
    rows = []
    for t in many["timeline_rows"]:
        cells = [("✓" if p else "·") for p in t["present"]]
        rows.append([t["gene"], t["label"]] + cells +
                    [t["first_seen"], "是" if t["fixed"] else "否"])
    parts.append(_table(head, rows))

    parts.append("<h2>各样本变异负荷与免疫命中</h2>")
    rows = []
    for name, r in many["per_sample"]:
        n_epi = sum(1 for x in r["impact_rows"] if x["in_epitope"])
        n_key = sum(1 for x in r["impact_rows"] if x["in_key_site"])
        rows.append([name, r["variant_summary"]["snp"], r["variant_summary"]["del"] + r["variant_summary"]["ins"],
                     len(r["impact_rows"]), n_epi, n_key,
                     ", ".join(m["lineage"] for m in r["lineage_matches"]) or "-"])
    parts.append(_table(["样本", "SNP", "indel事件", "蛋白突变", "落表位", "关键位点", "谱系签名"], rows))
    parts.append("<p class='muted'>✓=该突变存在于该时间点。固定=在最后一个时间点仍存在。"
                 "本报告为科研辅助输出，表位/结构影响为规则推断，非实验结论。</p>")
    parts.append("</body></html>")
    return "".join(parts)


def build_forecast_html(forecast: dict, meta: dict | None = None) -> str:
    css = """
    body{font-family:-apple-system,'Segoe UI',Arial,sans-serif;margin:28px;color:#1d2733}
    h1{font-size:22px} h2{font-size:17px;margin-top:26px;border-left:4px solid #d97706;padding-left:8px}
    table{border-collapse:collapse;width:100%;font-size:12.5px;margin:8px 0}
    th,td{border:1px solid #d4dbe3;padding:5px 8px;text-align:left;vertical-align:top}
    th{background:#fef3e2} .hi{color:#b91c1c;font-weight:700} .mid{color:#b45309}
    .muted{color:#64748b;font-size:12px}
    """
    parts = [f"<html><head><meta charset='utf-8'><style>{css}</style></head><body>"]
    parts.append("<h1>VirInsight 未来突变方向监测预测报告</h1>")
    parts.append(f"<p class='muted'>生成时间 {datetime.datetime.utcnow().isoformat(timespec='seconds')}Z · "
                 f"蛋白 {forecast['protein']} · 实验结构来源："
                 f"{forecast['structure']['pdb'] if forecast['structure'] else '无'} · "
                 f"实时频率：{'已接入' if forecast['live_ok'] else '不可用（使用知识快照）'}</p>")
    rows = []
    for x in forecast["rows"]:
        cls = "hi" if x["confidence"] == "高" else "mid" if x["confidence"] == "中" else ""
        rows.append([f"<span class='{cls}'>{x['mutation']}</span>", x["score"],
                     f"<span class='{cls}'>{x['confidence']}</span>",
                     x["direction"], x.get("rel_sasa"), "是" if x["interface"] else "",
                     x["blosum62"], ", ".join(x["seen_lineages"]) or "前瞻位点",
                     "；".join(x["reasons"][:3])])
    parts.append("<h2>重点监测突变优先级（多因子打分）</h2>")
    parts.append(_table(["候选突变", "评分", "置信度", "预测方向", "relSASA", "受体界面",
                         "BLOSUM62", "已见于", "判定依据（实验结构/知识快照/推断）"], rows))
    parts.append("<h2>因子构成</h2>")
    parts.append("<p class='muted'>F1 免疫选择压力(0-26) · F2 结构可及性SASA(0-12) · "
                 "F3 受体界面(0-16) · F4 生化可容忍度BLOSUM(0-12) · F5 流行率/逃逸先验(0-26) · "
                 "F6 用户时间序列轨迹(0-8)。</p>")
    frows = [[x["mutation"], x["factors"]["F1_immune"], x["factors"]["F2_accessibility"],
              x["factors"]["F3_interface"], x["factors"]["F4_biochem"],
              x["factors"]["F5_prior"], x["factors"]["F6_trajectory"]] for x in forecast["rows"]]
    parts.append(_table(["突变", "F1", "F2", "F3", "F4", "F5", "F6"], frows))
    parts.append(f"<p class='muted'>{html.escape(forecast['caveat'])}</p>")
    parts.append("<p class='muted'>免责声明：本预测为风险排序而非确定性预言，病毒突变存在强随机性。"
                 "结构因子来自实验PDB，逃逸判定含知识库与启发式成分，不构成临床结论。</p>")
    parts.append("</body></html>")
    return "".join(parts)
