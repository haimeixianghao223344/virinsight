"""基于实验 PDB 的三维结构环境评估。

对免疫重点蛋白（默认 SARS-CoV-2 S RBD / 6M0J）的每个残基提供真实三维量：
- 相对溶剂可及面积 relSASA（Bio.PDB Shrake-Rupley，Tien 2013 最大面积归一化）
- 配体/受体界面接触（重原子 ≤ 4.5Å；6M0J 中为 RBD–ACE2）
- 8Å 结构邻域残基
并对突变给出结构后果分类与"伪ΔΔG"启发式分（明确标注非 FoldX/Rosetta）。

序列映射：用 PDB 链自身的残基编号 + 序列比对，建立 S 蛋白全长位点 -> PDB 残基的映射。
"""
from __future__ import annotations
from pathlib import Path
import urllib.request
import warnings

import numpy as np

warnings.filterwarnings("ignore")

from Bio.PDB import PDBParser, ShrakeRupley
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1

PDB_DIR = Path(__file__).resolve().parent.parent / "data" / "pdb"

# 经验最大 SASA (Tien et al. 2013)
MAXASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "Q": 225,
          "E": 223, "G": 104, "H": 224, "I": 197, "L": 201, "K": 236,
          "M": 224, "F": 240, "P": 159, "S": 155, "T": 172, "W": 285,
          "Y": 263, "V": 174}

# 标准 BLOSUM62（仅放原型用到的；来自 NCBI 公开矩阵）
BLOSUM62 = {
 ("A","A"):4,("A","R"):-1,("A","N"):-2,("A","D"):-2,("A","C"):0,("A","Q"):-1,
 ("A","E"):-1,("A","G"):0,("A","H"):-2,("A","I"):-1,("A","L"):-1,("A","K"):-1,
 ("A","M"):-1,("A","F"):-2,("A","P"):-1,("A","S"):1,("A","T"):0,("A","W"):-3,
 ("A","Y"):-2,("A","V"):0,
 ("R","R"):5,("R","N"):0,("R","D"):-2,("R","C"):-3,("R","Q"):1,("R","E"):0,
 ("R","G"):-2,("R","H"):0,("R","I"):-3,("R","L"):-2,("R","K"):2,("R","M"):-1,
 ("R","F"):-3,("R","P"):-2,("R","S"):-1,("R","T"):-1,("R","W"):-3,("R","Y"):-2,
 ("R","V"):-3,
 ("N","N"):6,("N","D"):1,("N","C"):-3,("N","Q"):0,("N","E"):0,("N","G"):0,
 ("N","H"):1,("N","I"):-3,("N","L"):-3,("N","K"):0,("N","M"):-2,("N","F"):-3,
 ("N","P"):-2,("N","S"):1,("N","T"):0,("N","W"):-4,("N","Y"):-2,("N","V"):-3,
 ("D","D"):6,("D","C"):-3,("D","Q"):0,("D","E"):2,("D","G"):-1,("D","H"):-1,
 ("D","I"):-3,("D","L"):-4,("D","K"):-1,("D","M"):-3,("D","F"):-3,("D","P"):-1,
 ("D","S"):0,("D","T"):-1,("D","W"):-4,("D","Y"):-3,("D","V"):-3,
 ("C","C"):9,("C","Q"):-3,("C","E"):-4,("C","G"):-3,("C","H"):-3,("C","I"):-1,
 ("C","L"):-1,("C","K"):-3,("C","M"):-1,("C","F"):-2,("C","P"):-3,("C","S"):-1,
 ("C","T"):-1,("C","W"):-2,("C","Y"):-2,("C","V"):-1,
 ("Q","Q"):5,("Q","E"):2,("Q","G"):-2,("Q","H"):0,("Q","I"):-3,("Q","L"):-2,
 ("Q","K"):1,("Q","M"):0,("Q","F"):-3,("Q","P"):-1,("Q","S"):0,("Q","T"):-1,
 ("Q","W"):-2,("Q","Y"):-1,("Q","V"):-2,
 ("E","E"):5,("E","G"):-2,("E","H"):0,("E","I"):-3,("E","L"):-3,("E","K"):1,
 ("E","M"):-2,("E","F"):-3,("E","P"):-1,("E","S"):0,("E","T"):-1,("E","W"):-3,
 ("E","Y"):-2,("E","V"):-2,
 ("G","G"):6,("G","H"):-2,("G","I"):-4,("G","L"):-4,("G","K"):-2,("G","M"):-3,
 ("G","F"):-3,("G","P"):-2,("G","S"):0,("G","T"):-2,("G","W"):-2,("G","Y"):-3,
 ("G","V"):-3,
 ("H","H"):8,("H","I"):-3,("H","L"):-3,("H","K"):-1,("H","M"):-2,("H","F"):-1,
 ("H","P"):-2,("H","S"):-1,("H","T"):-2,("H","W"):-2,("H","Y"):2,("H","V"):-3,
 ("I","I"):4,("I","L"):2,("I","K"):-3,("I","M"):1,("I","F"):0,("I","P"):-3,
 ("I","S"):-2,("I","T"):-1,("I","W"):-3,("I","Y"):-1,("I","V"):3,
 ("L","L"):4,("L","K"):-2,("L","M"):2,("L","F"):0,("L","P"):-3,("L","S"):-2,
 ("L","T"):-1,("L","W"):-2,("L","Y"):-1,("L","V"):1,
 ("K","K"):5,("K","M"):-1,("K","F"):-3,("K","P"):-1,("K","S"):0,("K","T"):-1,
 ("K","W"):-3,("K","Y"):-2,("K","V"):-2,
 ("M","M"):5,("M","F"):0,("M","P"):-2,("M","S"):-1,("M","T"):-1,("M","W"):-1,
 ("M","Y"):-1,("M","V"):1,
 ("F","F"):6,("F","P"):-4,("F","S"):-2,("F","T"):-2,("F","W"):1,("F","Y"):3,
 ("F","V"):-1,
 ("P","P"):7,("P","S"):-1,("P","T"):-1,("P","W"):-4,("P","Y"):-3,("P","V"):-2,
 ("S","S"):4,("S","T"):1,("S","W"):-3,("S","Y"):-2,("S","V"):-2,
 ("T","T"):5,("T","W"):-2,("T","Y"):-2,("T","V"):0,
 ("W","W"):11,("W","Y"):2,("W","V"):-3,
 ("Y","Y"):7,("Y","V"):-1,
}


def blosum(a: str, b: str) -> int:
    if a == b:
        return BLOSUM62.get((a, a), 4)
    return BLOSUM62.get((a, b), BLOSUM62.get((b, a), -2))


def ensure_pdb(pdb_id: str) -> Path:
    PDB_DIR.mkdir(parents=True, exist_ok=True)
    path = PDB_DIR / f"{pdb_id.lower()}.pdb"
    if path.exists() and path.stat().st_size > 10000:
        return path
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "VirInsight/0.2"})
        path.write_bytes(urllib.request.urlopen(req, timeout=30).read())
    except Exception:  # noqa: BLE001
        if not path.exists():
            raise
    return path


# 结构注册表：病毒蛋白 -> 使用的PDB（含界面受体链与重点链）
STRUCTURE_REGISTRY = {
    ("sarscov2", "S"): {"pdb": "6m0j", "focus_chain": "E", "partner_chain": "A",
                        "partner_name": "ACE2", "covered_range": (333, 526),
                        "region": "RBD（受体结合结构域）"},
}


class StructureEnvironment:
    def __init__(self, pdb_id, focus_chain, partner_chain=None, partner_name=""):
        path = ensure_pdb(pdb_id)
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, str(path))
        ShrakeRupley().compute(structure, level="R")
        model = next(structure.get_models())
        self.pdb_id = pdb_id
        self.chain = model[focus_chain]
        self.partner = model[partner_chain] if partner_chain and partner_chain in model else None
        self.partner_name = partner_name
        self.residues = [r for r in self.chain if r.id[0] == " " and is_aa(r)]
        self.res_by_num = {r.id[1]: r for r in self.residues}
        self._seq = "".join(protein_letters_3to1[r.get_resname()] for r in self.residues)
        self.first_num = self.residues[0].id[1]
        self.interface = self._compute_interface()
        self._coords, self._resnums = self._heavy_table()

    def _heavy_table(self):
        coords, nums = [], []
        for r in self.residues:
            for a in r:
                if a.element != "H":
                    coords.append(a.coord)
                    nums.append(r.id[1])
        return np.array(coords), np.array(nums)

    def _compute_interface(self, cutoff=4.5):
        if self.partner is None:
            return set()
        pcoords = np.array([a.coord for r in self.partner if r.id[0] == " "
                            for a in r if a.element != "H"])
        interface = set()
        for r in self.residues:
            rc = np.array([a.coord for a in r if a.element != "H"])
            if len(rc) == 0:
                continue
            d = np.linalg.norm(pcoords[:, None, :] - rc[None, :, :], axis=2)
            if float(d.min()) <= cutoff:
                interface.add(r.id[1])
        return interface

    @property
    def sequence(self) -> str:
        """PDB链序列（顺序按残基号）。"""
        ordered = sorted(self.residues, key=lambda r: r.id[1])
        return "".join(protein_letters_3to1[r.get_resname()] for r in ordered)

    def residue_sasa(self, full_pos: int):
        """full_pos 为蛋白全长位点。PDB 编号 == S 蛋白编号（6M0J E 链 333-526）。"""
        r = self.res_by_num.get(full_pos)
        if r is None:
            return None
        aa = protein_letters_3to1[r.get_resname()]
        rel = float(r.sasa) / MAXASA.get(aa, 200)
        return dict(resnum=full_pos, aa=aa, abs_sasa=round(float(r.sasa), 1),
                    rel_sasa=round(rel, 3),
                    burial=("buried" if rel < 0.09 else "partial" if rel < 0.35 else "surface"),
                    interface=full_pos in self.interface)

    def neighborhood(self, full_pos: int, cutoff=8.0):
        r = self.res_by_num.get(full_pos)
        if r is None:
            return []
        rc = np.array([a.coord for a in r if a.element != "H"])
        if len(rc) == 0:
            return []
        d = np.linalg.norm(self._coords[:, None, :] - rc[None, :, :], axis=2).min(axis=1)
        near = sorted({int(self._resnums[i]) for i in np.where(d < cutoff)[0]
                       if int(self._resnums[i]) != full_pos})
        return near


def assess_mutation(env: StructureEnvironment, full_pos: int, ref_aa: str, alt_aa: str) -> dict:
    """返回该突变在真实结构环境下的后果评估；位点不在覆盖结构中则返回 covered=False。"""
    info = env.residue_sasa(full_pos)
    if info is None or info["aa"] != ref_aa:
        return dict(covered=False)
    near = env.neighborhood(full_pos)
    bscore = blosum(ref_aa, alt_aa)
    rel = info["rel_sasa"]
    interface = info["interface"]

    # 伪 ΔΔG：埋藏程度越高、BLOSUM越差，失稳风险越大（启发式，非物理计算）
    # 0(轻微) ~ 3(强烈失稳风险)
    destabilize = 0.0
    if bscore <= -3:
        destabilize += 2.2 * (1 - rel) + 0.6
    elif bscore == -2:
        destabilize += 1.4 * (1 - rel) + 0.4
    elif bscore == -1:
        destabilize += 0.7 * (1 - rel) + 0.2
    # 脯氨酸/甘氨酸进出核心的构象惩罚
    if (ref_aa in "PG") != (alt_aa in "PG") and rel < 0.25:
        destabilize += 0.8
    pseudo_ddg_band = min(3, int(round(destabilize)))

    reasons = []
    if interface:
        reasons.append(f"位于{env.partner_name}结合界面（重原子≤4.5Å接触），替换可直接改变受体结合亲和力")
    if info["burial"] == "buried":
        reasons.append("位于蛋白疏水核心（relSASA<0.09），侧链替换主要影响折叠稳定性")
    elif info["burial"] == "partial":
        reasons.append("部分埋藏（0.09≤relSASA<0.35），同时影响局部堆积与表面性质")
    else:
        reasons.append("溶剂暴露表面（relSASA≥0.35），对折叠稳定性影响通常较小，更可能影响抗原性/相互作用")
    if bscore <= -2:
        reasons.append(f"生化性质差异较大（BLOSUM62={bscore}）")
    elif bscore >= 1:
        reasons.append(f"保守替换（BLOSUM62={bscore}），结构后果预计较温和")
    if (ref_aa in "KRH" and alt_aa in "DE") or (ref_aa in "DE" and alt_aa in "KRH"):
        reasons.append("电荷翻转，可能改变界面盐桥与局部静电势")
    if (ref_aa in "PG") != (alt_aa in "PG"):
        reasons.append("脯氨酸/甘氨酸改变，可能影响局部主链构象")

    return dict(covered=True, ref_aa=ref_aa, alt_aa=alt_aa, pos=full_pos,
                rel_sasa=info["rel_sasa"], burial=info["burial"],
                interface=interface, blosum62=bscore,
                pseudo_ddg_band=pseudo_ddg_band,
                neighborhood=near[:14],
                partner_name=env.partner_name,
                reasons=reasons,
                evidence="experimental_structure+heuristic",
                caveat="伪ΔΔG为BLOSUM×埋藏度启发式，非FoldX/Rosetta物理自由能计算")
