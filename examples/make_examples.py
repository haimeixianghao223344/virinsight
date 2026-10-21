"""生成演示用 FASTA：在 Wuhan-Hu-1 参考基因组上按真实密码子位置引入各关注株
标志性氨基酸突变。这些是【合成演示序列】，用于在无法联网批量下载真实序列时
跑通流程；实际研究请上传真实测序基因组。

用法: python make_examples.py -> examples/demo_sarscov2_timeseries.fasta
"""
from __future__ import annotations
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio.Seq import Seq  # noqa: E402
from virinsight.core.io import load_reference, get_virus  # noqa: E402

# (基因, 氨基酸位, 目标氨基酸或'del') —— 选取知识库覆盖的标志性突变
SAMPLES = [
    ("2020-03_早期株(D614G)", [
        ("S", 614, "G"),
    ]),
    ("2021-06_Delta样(L452R+T478K+D614G+P681R)", [
        ("S", 452, "R"), ("S", 478, "K"), ("S", 614, "G"),
        ("S", 681, "R"), ("S", 19, "R"), ("S", 950, "S" if False else None),
        ("N", 203, "K"), ("N", 204, "R"),
    ]),
    ("2022-02_Omicron_BA1样(RBD多点逃逸)", [
        ("S", 67, "V"), ("S", 339, "D"), ("S", 371, "L"), ("S", 375, "F"),
        ("S", 417, "N"), ("S", 440, "K"), ("S", 477, "N"), ("S", 478, "K"),
        ("S", 496, "S"), ("S", 498, "R"), ("S", 501, "Y"), ("S", 505, "H"),
        ("S", 614, "G"), ("S", 655, "Y"), ("S", 681, "H"),
        ("N", 203, "K"), ("N", 204, "R"),
    ]),
    ("2023-04_XBB样(F486P+K444T等强逃逸)", [
        ("S", 180, "V"), ("S", 444, "T"), ("S", 456, "L"), ("S", 460, "K"),
        ("S", 478, "K"), ("S", 486, "P"), ("S", 501, "Y"),
        ("S", 614, "G"), ("S", 681, "H"),
        ("N", 203, "K"), ("N", 204, "R"),
    ]),
]

CODONS = {
    "A": ["GCT"], "R": ["AGA", "CGT"], "N": ["AAT"], "D": ["GAT"],
    "C": ["TGT"], "Q": ["CAA"], "E": ["GAA"], "G": ["GGT"], "H": ["CAT"],
    "I": ["ATT"], "L": ["TTA", "CTT"], "K": ["AAA"], "M": ["ATG"],
    "F": ["TTT"], "P": ["CCT"], "S": ["TCT", "AGT"], "T": ["ACT"],
    "W": ["TGG"], "Y": ["TAT"], "V": ["GTT"],
}


def best_codon(ref_codon: str, target_aa: str, table: int = 1) -> str:
    """选与原密码子差异最小的同义目标密码子。"""
    target_aa = target_aa.upper()
    cands = CODONS[target_aa]
    return min(cands, key=lambda c: sum(1 for a, b in zip(ref_codon, c) if a != b))


def main():
    import json
    ref = load_reference("sarscov2")
    v = get_virus("sarscov2")
    genes = {g["name"]: g for g in json.loads((Path(__file__).resolve().parent.parent /
             "virinsight" / "data" / "kb" / v["genes_kb"]).read_text(encoding="utf-8"))}
    out_lines = []
    log = []
    for name, muts in SAMPLES:
        s = list(ref)
        for gene, aapos, alt in muts:
            if alt is None:
                continue
            g = genes[gene]
            gstart0 = g["nt_start"] - 1 + (g.get("phase", 1) - 1)
            cstart = gstart0 + (aapos - 1) * 3
            ref_codon = ref[cstart:cstart + 3]
            ref_aa = str(Seq(ref_codon).translate())
            new_codon = best_codon(ref_codon, alt)
            new_aa = str(Seq(new_codon).translate())
            assert new_aa == alt, (gene, aapos, ref_codon, new_codon, new_aa)
            for i, b in enumerate(new_codon):
                s[cstart + i] = b
            log.append(f"{name}\t{gene}:{ref_aa}{aapos}{alt}\t{ref_codon}>{new_codon}\tnt{cstart+1}-{cstart+3}")
        out_lines.append(f">{name} | SYNTHETIC DEMO: derived from NC_045512.2 with engineered signature mutations")
        seq = "".join(s)
        for i in range(0, len(seq), 70):
            out_lines.append(seq[i:i + 70])
    outpath = Path(__file__).resolve().parent / "demo_sarscov2_timeseries.fasta"
    outpath.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    Path(__file__).resolve().parent.joinpath("demo_mutation_log.tsv").write_text(
        "sample\tprotein_change\tcodon\tnt_range\n" + "\n".join(log) + "\n", encoding="utf-8")
    print(f"written {outpath} ({len(SAMPLES)} samples)")
    print(f"mutations: {len(log)}")


if __name__ == "__main__":
    main()
