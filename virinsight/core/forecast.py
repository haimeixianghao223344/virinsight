"""未来突变方向预测（可解释的重点监测位点排序，而非确定性预言）。

多因子加权：
  F1 免疫选择压力  —— 表位/关键位点(知识库) + 邻近逃逸突变密度
  F2 结构可及性    —— 实验PDB relSASA（表面才可能被抗体接触）
  F3 受体界面      —— 6M0J RBD–ACE2 接触（可改变结合/感染性）
  F4 生化可容忍度  —— BLOSUM62（越高=越可能被结构允许而保留）
  F5 流行率先验    —— 已在关注株观察到（知识快照）；可选实时频率
  F6 时间轨迹      —— 用户多时期样本中晚期新出现/累积位点
输出 0-100 监测优先级、置信度分级、分类方向（抗体逃逸/受体结合/T细胞），
并对每条候选给人类可读依据。所有因子透明列出。
"""
from __future__ import annotations
from pathlib import Path
import json

from .io import KB_DIR
from .structure_impact import (StructureEnvironment, assess_mutation,
                               STRUCTURE_REGISTRY, blosum)
from .immuno import residue_class
from . import frequency as freq_mod


def load_prior(name: str = "sarscov2_prior.json") -> dict:
    return json.loads((KB_DIR / name).read_text(encoding="utf-8"))


def _epitope_and_key(pos: int, epi_kb: dict, protein: str = "S"):
    """返回(是否在表位, 是否关键位点, 描述列表)。"""
    pkb = epi_kb.get("proteins", {}).get(protein, {})
    desc = []
    in_epi = False
    for bucket, lab in [("bcell_neutralizing", "中和/B细胞表位"),
                        ("tcell_mhc1", "MHC-I T细胞表位"),
                        ("tcell_mhc2", "MHC-II T细胞表位")]:
        for ep in pkb.get(bucket, []):
            if ep["aa_start"] <= pos <= ep["aa_end"]:
                in_epi = True
                desc.append(f"{lab}:{ep['name']}")
    is_key = False
    for kr in pkb.get("key_residues", []):
        if kr["aa"] == pos:
            is_key = True
            desc.append(f"关键位点 {kr['label']}({kr['aspect']})")
    return in_epi, is_key, desc


def forecast_sites(virus_id: str, protein: str, epi_kb: dict,
                   timeseries=None, fetch_live: bool = False,
                   top_n: int = 20) -> dict:
    """生成监测位点列表。
    timeseries: 可选 analyze_many 结果（用于F6轨迹信号）。
    """
    reg = STRUCTURE_REGISTRY.get((virus_id, protein))
    env = None
    if reg:
        try:
            env = StructureEnvironment(reg["pdb"], reg["focus_chain"],
                                       reg["partner_chain"], reg["partner_name"])
        except Exception:  # noqa: BLE001
            env = None

    prior = load_prior()
    pos_records = prior["positions"]
    # 已在用户样本晚期出现的位点（F6）
    trajectory = {}
    if timeseries is not None:
        names = timeseries["sample_names"]
        for t in timeseries["timeline_rows"]:
            if t["gene"] == protein and any(t["present"]):
                first = next((i for i, p in enumerate(t["present"]) if p), 0)
                trajectory[t["aa_pos"]] = dict(
                    first_idx=first, late=first >= max(0, len(names) - 2),
                    label=t["label"])

    # 实时频率（best-effort）
    live = {}
    if fetch_live and env is not None:
        mut_codes = [f"{protein}:{p['pos']}{p['candidate_alts'][0]}"
                     for p in pos_records[:12]]
        try:
            live = freq_mod.fetch_live_prior(mut_codes)
        except Exception:  # noqa: BLE001
            live = {}

    rows = []
    for rec in pos_records:
        pos = rec["pos"]
        in_epi, is_key, edesc = _epitope_and_key(pos, epi_kb, protein)
        seen = rec.get("seen_lineages", [])
        seen_strength = 3 if len(seen) >= 3 else 2 if seen else 0

        # 参考氨基酸：优先取实验结构中的真实残基；无结构时退回知识库标签
        if env is not None:
            _info = env.residue_sasa(pos)
            ref_aa = _info["aa"] if _info else (_ref_aa(pos, prior, epi_kb))
        else:
            ref_aa = _ref_aa(pos, prior, epi_kb)
        # 每个候选替换分别评估，取该位点最强（最值得监测）的一个
        best = None
        for alt in rec["candidate_alts"]:
            if alt == ref_aa:
                continue  # 同义/未改变，跳过
            factors = {}
            # F1 免疫压力 0-1
            f1 = 0.0
            if is_key:
                f1 = max(f1, 1.0)
            elif in_epi:
                f1 = max(f1, 0.75)
            factors["F1_immune"] = round(f1, 2)

            # F2/F3 结构
            rel_sasa = burial = None
            interface = False
            bl = blosum(ref_aa, alt)
            f2 = f3 = 0.0
            struct_note = []
            if env is not None:
                a = assess_mutation(env, pos, ref_aa, alt)
                if a["covered"]:
                    rel_sasa, burial = a["rel_sasa"], a["burial"]
                    interface = a["interface"]
                    f2 = min(1.0, rel_sasa) if rel_sasa is not None else 0.3
                    f3 = 1.0 if interface else 0.0
                    bl = a["blosum62"]
                    struct_note = a["reasons"]
            else:
                f2 = 0.4  # 无结构时中性
            factors["F2_accessibility"] = round(f2, 2)
            factors["F3_interface"] = round(f3, 2)

            # F4 生化可容忍度（BLOSUM -4..11 -> 0..1，越可容忍越高）
            f4 = max(0.0, min(1.0, (bl + 4) / 12))
            factors["F4_biochem"] = round(f4, 2)

            # F5 先验
            f5 = min(1.0, (seen_strength / 3.0) * 0.6 + rec.get("escape_strength", 0) / 3 * 0.4)
            factors["F5_prior"] = round(f5, 2)

            # F6 轨迹
            f6 = 0.6 if trajectory.get(pos, {}).get("late") else 0.0
            factors["F6_trajectory"] = round(f6, 2)

            live_note = ""
            live_key = f"{protein}:{pos}{alt}"
            if live_key in live:
                f5 = max(f5, 0.8)
                factors["F5_prior"] = round(f5, 2)
                live_note = f"实时频率={live[live_key].get('prevalence')}"

            # 权重：免疫与逃逸先验主导，结构可及性/界面强化，可容忍度作为可行性约束
            score = (26 * f1 + 12 * f2 + 16 * f3 + 12 * f4 + 26 * f5 + 8 * f6)
            score = round(min(100, score), 1)

            direction = []
            if in_epi or is_key or rec.get("escape_strength", 0) >= 2:
                direction.append("抗体/T细胞逃逸")
            if interface:
                direction.append("受体结合/感染性")
            if not direction:
                direction.append("一般适应性")

            reasons = list(edesc)
            reasons.append(f"逃逸先验强度{rec.get('escape_strength','?')}/3"
                           + (f"，已见于{','.join(seen)}" if seen else "（关注株尚未广泛固定，属前瞻位点）"))
            if struct_note:
                reasons.extend(struct_note[:2])
            if pos in trajectory:
                reasons.append(f"用户时间序列中{trajectory[pos]['label']}于晚期出现")
            if live_note:
                reasons.append(live_note + f"（{freq_mod.UA['User-Agent'].split('/')[0]}实时）")

            conf = ("高" if (is_key and seen_strength >= 2) else
                    "中" if (in_epi or seen_strength >= 2 or interface) else "低")

            cand = dict(pos=pos, alt=alt, mutation=f"{ref_aa}{pos}{alt}",
                        score=score, confidence=conf,
                        direction="、".join(direction), factors=factors,
                        rel_sasa=rel_sasa, burial=burial, interface=interface,
                        blosum62=bl, seen_lineages=seen,
                        reasons=reasons,
                        ref_aa=ref_aa,
                        note=rec.get("note", ""))
            if best is None or cand["score"] > best["score"]:
                best = cand
        if best is not None:
            rows.append(best)

    rows.sort(key=lambda r: r["score"], reverse=True)
    return dict(protein=protein, structure=reg,
                structure_ok=env is not None,
                live_ok=bool(live),
                rows=rows[:top_n], all_rows=rows,
                caveat=("监测优先级为多因子统计/结构推断，非确定性预测；"
                        "结构因子来自实验PDB(6M0J RBD-ACE2)，流行率为知识快照"
                        + ("+实时API" if live else "（实时频率源当前不可用，已降级）")
                        + "；伪ΔΔG非FoldX/Rosetta物理计算。"))


def _ref_aa(pos: int, prior: dict, epi_kb: dict) -> str:
    """从知识库关键位点标签推断参考氨基酸；缺省返回空。"""
    pkb = epi_kb.get("proteins", {}).get("S", {})
    for kr in pkb.get("key_residues", []):
        if kr["aa"] == pos and kr.get("label"):
            return kr["label"][0]
    return "X"
