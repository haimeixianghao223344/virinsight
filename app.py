"""VirInsight — 病毒基因组免疫与突变分析平台 (Streamlit 原型 v0.1)
启动: streamlit run app.py
"""
from __future__ import annotations
import io
import warnings

warnings.filterwarnings("ignore")
import pandas as pd
import plotly.express as px
import streamlit as st

from virinsight.core.io import load_reference, get_virus, load_viruses, load_fasta, parse_fasta_text
from virinsight.core.pipeline import identify_virus, analyze_sample, analyze_many
from virinsight.report import build_single_html, build_timeseries_html, build_forecast_html
from virinsight.core.forecast import forecast_sites
from virinsight.core.io import load_json

st.set_page_config(page_title="VirInsight 病毒免疫与突变分析", layout="wide")

VIRUSES = {v["id"]: v for v in load_viruses()}
EXAMPLE_FASTA = "examples/demo_sarscov2_timeseries.fasta"


def _embed_3dmol(pdb_id: str, chain: str, highlight: list[int]):
    """内联3Dmol.js CDN；无网络环境下组件不渲染（表格仍可用）。"""
    import streamlit.components.v1 as components
    sel_hi = ",".join(str(p) for p in highlight)
    js = """
    <div id="vview" style="height:430px;width:100%"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
      try {
        var el = document.getElementById('vview');
        var viewer = $3Dmol.createViewer(el, {backgroundColor:'white'});
        $3Dmol.download('pdb:__PDB__', viewer, {}, function(){
          viewer.setStyle({chain:'__CHAIN__'}, {cartoon:{color:'#7aa7e8'}});
          viewer.setStyle({chain:'A'}, {cartoon:{color:'#dddddd',opacity:0.45}});
          viewer.setStyle({chain:'__CHAIN__',resi:[__HILITE__]},
            {cartoon:{color:'orange'},stick:{colorscheme:'yellowCarbon'}});
          viewer.zoomTo({chain:'__CHAIN__'});
          viewer.render();
        });
      } catch(e){ document.getElementById('vview').innerHTML='3D 组件加载失败（可能离线），请参考表格。'; }
    </script>"""
    js = js.replace("__PDB__", pdb_id).replace("__CHAIN__", chain).replace("__HILITE__", sel_hi)
    components.html(js, height=450)


def infer_date_key(rec_id: str, desc: str, idx: int):
    import re
    m = re.search(r"(20\d{2})[-_/]?(\d{2})?", desc)
    if m:
        return m.group(1) + "-" + (m.group(2) or "01")
    return f"样本{idx+1}"


def load_input(uploaded_file, pasted: str, use_example: bool):
    if use_example:
        recs = load_fasta(EXAMPLE_FASTA)
    elif uploaded_file is not None:
        recs = load_fasta(uploaded_file)
    elif pasted.strip():
        recs = parse_fasta_text(pasted)
    else:
        return []
    return [(rid.split()[0] if rid else f"seq{i}", desc, seq) for i, (rid, desc, seq) in enumerate(recs)]


# ============ 侧栏 ============
st.sidebar.title("VirInsight")
st.sidebar.caption("病毒基因组免疫与突变分析原型 v0.1 · 科研辅助")
tab = st.sidebar.radio("分析模式", ["① 单基因组免疫画像", "② 多时期突变追踪", "③ 蛋白结构影响(重点蛋白)", "④ 未来突变方向预测"])

st.sidebar.markdown("---")
auto_id = st.sidebar.checkbox("自动识别病毒种类", value=True)
vid_manual = st.sidebar.selectbox("或手动指定参考库", list(VIRUSES.keys()),
                                  format_func=lambda k: VIRUSES[k]["name"], disabled=auto_id)
fetch_online = st.sidebar.checkbox("在线增强注释 (UniProt/RCSB，无网自动降级)", value=False)
st.sidebar.markdown("---")
up = st.sidebar.file_uploader("上传 FASTA（1条=单样本，多条=时间序列）", type=["fa", "fasta", "fna", "txt"])
pasted = st.sidebar.text_area("或粘贴 FASTA 序列", height=120)
use_example = st.sidebar.checkbox("使用内置示例（SARS-CoV-2 四时期合成演示序列）", value=True)

recs = load_input(up, pasted, use_example and up is None and not pasted.strip())
if not recs:
    st.info("请在左侧上传/粘贴 FASTA，或勾选内置示例。")
    st.stop()

# 识别病毒
first_seq = recs[0][2]
if auto_id:
    with st.spinner("正在识别病毒种类…"):
        cands = identify_virus(first_seq)
    top = cands[0]
    if top[3] < 0.15:
        st.warning(f"未在内置参考库找到近缘病毒（最高相似度仅 {top[3]:.1%}）。"
                   "可手动选择参考库运行；非内置病毒将只做通用变异与结构域分析。")
    vid = top[0]
    with st.sidebar.expander(f"识别结果：{VIRUSES[vid]['name']}（相似度 {top[3]:.1%}）"):
        st.caption("点击展开各候选")
        for cid, cname, jac, frac, _ in cands:
            st.write(f"- {cname}: k-mer 覆盖 {frac:.1%}")
else:
    vid = vid_manual
virus = VIRUSES[vid]
ref = load_reference(vid)
st.sidebar.success(f"参考株：{virus['name']}\n\n{virus['ref_name']}")

online_cache = st.session_state.setdefault("online_cache", {})

# ============ 模式① 单基因组 ============
if tab.startswith("①"):
    st.title("单基因组免疫画像")
    if len(recs) > 1:
        idx = st.selectbox("选择要分析的样本", range(len(recs)),
                           format_func=lambda i: recs[i][0])
    else:
        idx = 0
    rid, desc, seq = recs[idx]
    st.caption(f"样本 `{rid}` · 长度 {len(seq)} nt · {desc[:120]}")
    with st.spinner("比对、翻译并注释免疫位点…"):
        res = analyze_sample(ref, seq, virus, fetch_online=fetch_online, online_cache=online_cache)
    vs = res["variant_summary"]
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("SNP", vs["snp"])
    c2.metric("缺失事件", f"{vs['del']} ({vs['del_nt']}nt)")
    c3.metric("插入事件", f"{vs['ins']} ({vs['ins_nt']}nt)")
    c4.metric("命中免疫/关键位点", sum(1 for r in res["impact_rows"] if r["in_epitope"] or r["in_key_site"]))

    st.subheader("免疫主作用基因排序")
    st.caption("综合免疫相关度（知识库★）与该样本在表位/关键位点上的实际突变；类别含抗原/免疫拮抗/复制酶等")
    gr = pd.DataFrame([{"基因": g["name"], "类别": g["category"],
                        "免疫相关度": "★" * g["immune_relevance"],
                        "错义突变": g["n_missense"], "蛋白indel": g["n_indel"],
                        "命中免疫位点": g["n_immuno_hits"],
                        "免疫作用": g["immuno_note"]}
                       for g in res["ranked_genes"]
                       if g["n_missense"] + g["n_indel"] + g["n_immuno_hits"] > 0])
    st.dataframe(gr, use_container_width=True, hide_index=True)

    if res["lineage_matches"]:
        st.subheader("关注株/谱系签名匹配")
        st.dataframe(pd.DataFrame([{"谱系": m["lineage"], "命中": f"{m['hit']}/{m['total']}",
                                    "比例": m["fraction"], "命中位点": ", ".join(m["matched"])}
                                   for m in res["lineage_matches"]]),
                     use_container_width=True, hide_index=True)

    st.subheader("蛋白突变 → 结构与免疫影响")
    ir = pd.DataFrame([{"基因": r["gene"], "突变": r["label"],
                        "类型": r["type"], "结构严重性(0-3)": r["severity"],
                        "关键位点": "是" if r["in_key_site"] else "",
                        "位于表位": "是" if r["in_epitope"] else "",
                        "影响判定": r["immune_tags"]} for r in res["impact_rows"]])
    st.dataframe(ir, use_container_width=True, hide_index=True)

    with st.expander("核酸变异（相对参考株）"):
        st.dataframe(pd.DataFrame([{"参考位置": v["ref_pos"], "类型": v["type"],
                                    "参考碱基": v["ref_nt"], "变异碱基": v["alt_nt"]}
                                   for v in res["variants"]]),
                     use_container_width=True, hide_index=True, height=260)

    html_bytes = build_single_html(dict(sample=rid, length=len(seq)), res, virus).encode("utf-8")
    st.download_button("下载 HTML 报告", html_bytes, f"virinsight_single_{rid}.html", "text/html")

# ============ 模式② 多时期 ============
elif tab.startswith("②"):
    st.title("多时期基因组突变追踪")
    if len(recs) < 2:
        st.warning("该模式需要多条序列（按 FASTA 中的顺序视为时间先后）。当前只有 1 条。")
        st.stop()
    # 尝试按描述中的日期排序，允许用户关闭
    sort_by_date = st.checkbox("按序列名中的日期自动排序", value=True)
    indexed = list(enumerate(recs))
    if sort_by_date:
        indexed.sort(key=lambda t: infer_date_key(recs[t[0]][0], recs[t[0]][1], t[0]))
    samples = [(recs[i][0], recs[i][2]) for i, _ in indexed]
    st.caption("时间顺序：" + " → ".join(n for n, _ in samples))
    with st.spinner("逐样本比对并构建突变时间线…"):
        many = analyze_many(ref, samples, virus)

    # 热图
    mat = {t["label"]: t["present"] for t in many["timeline_rows"]}
    genes = [t["gene"] for t in many["timeline_rows"]]
    labels = list(mat.keys())
    if labels:
        dfh = pd.DataFrame(mat, index=many["sample_names"]).T
        dfh.insert(0, "基因", genes)
        st.subheader("突变 × 时间点 热图")
        fig = px.imshow(dfh.drop(columns=["基因"]).values,
                        x=many["sample_names"], y=labels,
                        color_continuous_scale=[[0, "#eef2f7"], [1, "#16a34a"]],
                        aspect="auto", labels=dict(color="存在"))
        fig.update_traces(showscale=False, xgap=2, ygap=2)
        fig.update_layout(height=max(300, 26 * len(labels)), font=dict(size=11))
        st.plotly_chart(fig, use_container_width=True)

    st.subheader("突变时间线（含首次出现 / 末期是否固定）")
    rows = [{"基因": t["gene"], "突变": t["label"], **{n: ("✓" if p else "") for n, p in
             zip(many["sample_names"], t["present"])},
             "首次出现": t["first_seen"], "末期固定": "是" if t["fixed"] else ""}
            for t in many["timeline_rows"]]
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    st.subheader("各样本突变负荷与谱系签名")
    load_rows = []
    for name, r in many["per_sample"]:
        load_rows.append({"样本": name, "SNP": r["variant_summary"]["snp"],
                          "indel事件": r["variant_summary"]["del"] + r["variant_summary"]["ins"],
                          "蛋白突变": len(r["impact_rows"]),
                          "落表位": sum(1 for x in r["impact_rows"] if x["in_epitope"]),
                          "关键位点": sum(1 for x in r["impact_rows"] if x["in_key_site"]),
                          "谱系签名": ", ".join(m["lineage"] for m in r["lineage_matches"]) or "-"})
    st.dataframe(pd.DataFrame(load_rows), use_container_width=True, hide_index=True)
    fig2 = px.bar(pd.DataFrame(load_rows), x="样本",
                  y=["SNP", "蛋白突变", "落表位", "关键位点"], barmode="group",
                  title="各时期突变负荷")
    st.plotly_chart(fig2, use_container_width=True)

    html_bytes = build_timeseries_html(dict(), many, virus).encode("utf-8")
    st.download_button("下载 HTML 时间线报告", html_bytes,
                       "virinsight_timeseries.html", "text/html")

# ============ 模式③ 结构影响 ============
elif tab.startswith("③"):
    st.title("蛋白结构影响（重点蛋白）")
    st.caption("把突变映射到功能结构域 / 受体界面 / 糖基化位点，并关联实验 PDB 结构。"
               "结构严重性为启发式；ΔΔG 计算（FoldX/Rosetta）为预留接口。")
    rid, desc, seq = recs[0] if len(recs) == 1 else (recs[0][0], recs[0][1], recs[0][2])
    res = analyze_sample(ref, seq, virus, fetch_online=fetch_online, online_cache=online_cache)

    st.subheader("已知实验结构（内置参考）")
    st.dataframe(pd.DataFrame(virus.get("pdbs", [])), use_container_width=True, hide_index=True)

    focus_options = [g["name"] for g in res["ranked_genes"]
                     if res["annotated"][g["name"]]["mutations"]]
    if not focus_options:
        st.info("该样本未检出蛋白突变。")
    focus = st.selectbox("选择重点蛋白", focus_options or ["S"])
    gdef = next(g for g in res["ranked_genes"] if g["name"] == focus)
    ann = res["annotated"][focus]

    st.markdown(f"**{focus}** — {gdef['function']}")
    if gdef.get("domains"):
        st.markdown("结构域/功能区：")
        st.dataframe(pd.DataFrame(gdef["domains"]), use_container_width=True, hide_index=True)

    st.markdown("该蛋白突变的结构/免疫影响：")
    rows = [r for r in res["impact_rows"] if r["gene"] == focus]

    # 三维结构环境增强（若该蛋白已登记实验PDB）
    from virinsight.core.structure_impact import StructureEnvironment, assess_mutation, STRUCTURE_REGISTRY
    reg3d = STRUCTURE_REGISTRY.get((vid, focus))
    env3d = None
    if reg3d:
        try:
            env3d = StructureEnvironment(reg3d["pdb"], reg3d["focus_chain"],
                                         reg3d["partner_chain"], reg3d["partner_name"])
            st.success(f"已加载实验结构 {reg3d['pdb']}（{reg3d['region']}–{reg3d['partner_name']}复合物），"
                       f"对 {reg3d['covered_range'][0]}-{reg3d['covered_range'][1]} 位残基提供真实 SASA/界面评估")
        except Exception:
            env3d = None

    table_rows = []
    for r in rows:
        entry = {"突变": r["label"], "类型": r["type"],
                 "结构严重性(0-3)": r["severity"],
                 "关键位点": "是" if r["in_key_site"] else "",
                 "影响": r["immune_tags"]}
        if env3d is not None and r["type"] in ("missense", "nonsense") and len(r["ref_aa"]) == 1 and len(r["alt_aa"]) == 1:
            a = assess_mutation(env3d, r["aa_pos"], r["ref_aa"], r["alt_aa"])
            if a["covered"]:
                entry.update({"relSASA": a["rel_sasa"], "埋藏度": a["burial"],
                              "受体界面": "是" if a["interface"] else "",
                              "伪ΔΔG档(0-3)": a["pseudo_ddg_band"], "BLOSUM62": a["blosum62"]})
        table_rows.append(entry)
    st.dataframe(pd.DataFrame(table_rows), use_container_width=True, hide_index=True)
    if env3d is not None:
        with st.expander("三维结构后果逐条解读（实验结构 + 启发式）"):
            for r in rows:
                if r["type"] in ("missense", "nonsense") and len(r["ref_aa"]) == 1 and len(r["alt_aa"]) == 1:
                    a = assess_mutation(env3d, r["aa_pos"], r["ref_aa"], r["alt_aa"])
                    if a["covered"]:
                        st.markdown(f"**{r['label']}**（relSASA={a['rel_sasa']}, {a['burial']}, "
                                    f"{'受体界面' if a['interface'] else '非界面'}, "
                                    f"伪ΔΔG档={a['pseudo_ddg_band']}, BLOSUM62={a['blosum62']}）")
                        for why in a["reasons"]:
                            st.write("- " + why)
            st.caption("伪ΔΔG为BLOSUM×埋藏度启发式，非FoldX/Rosetta物理自由能；界面由6M0J中重原子≤4.5Å判定。")

    if fetch_online and focus in res.get("online", {}):
        a = res["online"][focus]
        u = a["uniprot"]
        st.markdown(f"**在线注释** UniProt {u['accession']}（{u.get('source')}，{u['fetched_at']}）")
        if u.get("ok"):
            st.dataframe(pd.DataFrame(u["features"]), use_container_width=True, hide_index=True)
            st.caption("RCSB 检索到的实验结构: " + ", ".join(a["rcsb"].get("structures", [])[:15]))
        else:
            st.caption(f"在线获取失败（{u.get('error')}），已使用本地知识库。")
    elif fetch_online:
        acc = virus.get("uniprot_focus", {}).get(focus)
        if acc:
            with st.spinner(f"从 UniProt/RCSB 获取 {focus} ({acc}) …"):
                from virinsight.core.online import annotate_online
                online_cache[acc] = annotate_online(acc)
            st.rerun()

# ============ 模式④ 未来突变方向预测 ============
elif tab.startswith("④"):
    st.title("未来突变方向预测（重点监测位点排序）")
    st.warning("病毒突变具强随机性，本页输出的是**可解释的风险排序(watchlist)而非确定性预言**。"
               "结构因子来自实验PDB，流行率先验为知识快照（实时源不可达时降级），请结合专业判断。")

    if vid != "sarscov2":
        st.info("预测模块的结构(6M0J)与先验库当前内置在 SARS-CoV-2 S 蛋白；"
                "其他病毒可使用前三页的通用变异/免疫/在线注释。以下仍展示知识库可用部分。")

    epi = load_json(virus["epitopes_kb"])
    # 时间序列轨迹信号（若上传多条）
    ts_obj = None
    if len(recs) >= 2:
        with st.spinner("计算多时期轨迹信号…"):
            _samples = [(r[0], r[2]) for r in recs]
            ts_obj = analyze_many(ref, _samples, virus)
        st.caption("已纳入时间轨迹因子：" + " → ".join(ts_obj["sample_names"]))

    use_live = st.checkbox("尝试接入实时全球频率 API（outbreak.info，失败自动降级）", value=False)
    with st.spinner("整合实验结构 / 表位 / 先验，计算监测优先级…"):
        fc = forecast_sites(vid, "S", epi, timeseries=ts_obj, fetch_live=use_live, top_n=20)

    c1, c2 = st.columns(2)
    c1.metric("候选监测突变", len(fc["all_rows"]))
    c2.metric("实验结构", fc["structure"]["pdb"] if fc["structure"] else "无（降级）")
    if not fc["structure_ok"]:
        st.info("实验结构未能加载（可能无网络且本地无缓存PDB），F2/F3以中性先验替代。")
    if not fc["live_ok"]:
        st.caption("实时频率源当前不可用，F5 使用知识快照。")

    st.subheader("重点监测突变优先级")
    df = pd.DataFrame([{"候选突变": x["mutation"], "评分": x["score"], "置信度": x["confidence"],
                        "预测方向": x["direction"], "relSASA": x["rel_sasa"],
                        "埋藏度": x["burial"], "受体界面": "是" if x["interface"] else "",
                        "BLOSUM62": x["blosum62"],
                        "已见于谱系": ", ".join(x["seen_lineages"]) or "前瞻位点"}
                       for x in fc["rows"]])
    st.dataframe(df, use_container_width=True, hide_index=True)

    fig = px.bar([x for x in fc["rows"]][::-1], x="score", y="mutation",
                 orientation="h", color="confidence",
                 color_discrete_map={"高": "#dc2626", "中": "#d97706", "低": "#65a30d"},
                 title="监测优先级评分（Top 20）", labels={"score": "评分", "mutation": ""})
    fig.update_layout(height=560, font=dict(size=11))
    st.plotly_chart(fig, use_container_width=True)

    st.subheader("逐位点判定依据（透明因子）")
    for x in fc["rows"]:
        with st.expander(f"{x['mutation']} · 评分{x['score']} · {x['confidence']}置信 · {x['direction']}"):
            st.write("**方向**：" + x["direction"])
            st.write("**依据**：")
            for r in x["reasons"]:
                st.write("- " + r)
            fac = x["factors"]
            st.caption(f"F1免疫={fac['F1_immune']} F2可及性={fac['F2_accessibility']} "
                       f"F3界面={fac['F3_interface']} F4生化={fac['F4_biochem']} "
                       f"F5先验={fac['F5_prior']} F6轨迹={fac['F6_trajectory']}")
            if x.get("rel_sasa") is not None:
                st.caption(f"结构：relSASA={x['rel_sasa']}（{x['burial']}），"
                           f"{'位于受体界面' if x['interface'] else '非界面'}，BLOSUM62={x['blosum62']}")

    # 3D 视图（有 PDB 时；CDN 无网自动不渲染）
    if fc["structure_ok"]:
        with st.expander("在三维结构中查看监测位点（RBD–ACE2, 6M0J）"):
            st.caption("红=ACE2界面残基；橙=Top监测候选。需联网加载3Dmol.js，无网时请参考上方表格。")
            top_pos = [x["pos"] for x in fc["rows"][:12]]
            iface = sorted({p for x in fc["rows"] for p in []})  # 界面由结构内部决定
            _embed_3dmol("6m0j", "E", top_pos)

    html_bytes = build_forecast_html(fc).encode("utf-8")
    st.download_button("下载预测 HTML 报告", html_bytes,
                       "virinsight_forecast.html", "text/html")
    st.caption(fc["caveat"])

st.sidebar.markdown("---")
st.sidebar.caption("证据标注：★知识库精选(文献/数据库共识)｜启发式=规则推断。"
                   "表位逃逸≠实验验证；不构成临床结论。")
