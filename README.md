# VirInsight — 病毒基因组免疫与突变分析平台（原型 v0.2）

输入病毒基因组，自动识别病毒种类、注释基因、**定位与人体免疫相关的主要作用基因**（抗原 / 中和抗体表位 / T 细胞表位 / 免疫拮抗蛋白）；
输入**不同时期的多个基因组**，可追踪突变时间线、翻译为蛋白突变，并判定其落在**受体界面 / 中和表位 / T 细胞表位 / 功能结构域 / 糖基化位点**的结构与免疫影响。

> 定位：科研辅助原型。表位逃逸与"结构严重性"为**知识库规则 + 在线数据库推断**，不是 ΔΔG 计算或实验结果，不构成临床结论。

## 快速开始

```bash
pip install -r requirements.txt
streamlit run app.py
```

浏览器打开后，左侧默认勾选「使用内置示例（SARS-CoV-2 四时期合成演示序列）」即可一键体验。
也可以上传自己的 FASTA（1 条=单样本；多条按序列顺序/名称中的日期视为时间序列），或直接粘贴 FASTA。

- 勾选「在线增强注释」会实时访问 UniProt（活性位点/结构域/跨膜/糖基化）与 RCSB PDB（实验结构）；无网络时自动降级，本地流程不受影响。

## 三个分析页签

1. **单基因组免疫画像**：变异统计 → 免疫主作用基因排序 → 关注株/谱系签名匹配 → 蛋白突变的结构/免疫影响明细 → 核酸变异；可导出 HTML。
2. **多时期突变追踪**：突变 × 时间点热图、突变首次出现/末期固定、各时期突变负荷与谱系演变；可导出 HTML。
3. **蛋白结构影响（重点蛋白）**：突变映射到结构域/受体界面/糖基化位点，关联内置 PDB 与 UniProt 在线注释。

## 分析方法（核心链路）

| 步骤 | 方法 | 文件 |
|---|---|---|
| 病毒识别 | 参考库 21-mer k-mer 轮廓（Jaccard/样本覆盖度） | `core/identify.py` |
| 核酸变异 | edlib 全局比对（C，毫秒级）回溯 CIGAR；差异块用高 gap 罚分局部重比对，避免同聚区把相邻 SNP 误判为 indel | `core/variant.py` |
| 基因注释/翻译 | 知识库坐标提取 ORF，按相位翻译；氨基酸级比对生成 p.D614G 式突变 | `core/annotate.py` |
| 同义/非同义 | 密码子级统计 syn/nonsyn/stop（dN/dS 粗估，未位点标准化） | `core/annotate.py` |
| 免疫影响 | 突变映射到精选表位（MHC-I/II、中和/B细胞）、关键功能残基、结构域、furin 切位、N-糖基化获得/丢失；0–3 启发式结构严重性 | `core/immuno.py` |
| 谱系签名 | 与内置关注株标志突变集（Alpha/Delta/BA.1/XBB.1.5 等）做集合匹配 | 知识库 + `core/pipeline.py` |
| 在线增强 | UniProt REST、RCSB Search REST，失败优雅降级 | `core/online.py` |
| 三维结构影响 | PDB下载缓存、SASA、受体界面、邻域、伪ΔΔG | `core/structure_impact.py` |
| 未来突变预测 | F1-F6多因子监测位点排序 + 实时频率(best-effort) | `core/forecast.py`, `core/frequency.py` |

## 内置参考与知识库

- 参考基因组：SARS-CoV-2 Wuhan-Hu-1（NC_045512.2，NCBI 实拉）、HIV-1 HXB2（K03455.1，NCBI 实拉）。
- 知识库（`virinsight/data/kb/*.json`，可自行扩展，无需改代码）：
  - 基因功能、免疫相关度（★1–5）、免疫作用说明、结构域坐标、糖基化位点；
  - SARS-CoV-2：S 蛋白 RBD/RBM/NTD/切位，关键逃逸位点（K417/N439/K444/L452/T478/E484/F486/N501/Y505、D614、P681 等）、MHC-I/II 与中和表位；N/ORF3a；免疫拮抗 ORF6/ORF8；
  - HIV-1：Env gp120/gp41（V3/CD4bs/MPER/糖盾）、Gag p17/p24 经典 CTL 表位（SLYNTVATL、KK10、TW10）、Nef/Vif/Vpu 免疫拮抗、Pol 耐药位点。

**新增病毒**：在 `data/refs/` 放参考 FASTA，在 `data/kb/` 增加 `xxx_genes.json` / `xxx_epitopes.json`，并在 `viruses.json` 注册即可。
没有知识库的病毒仍可运行"通用变异 + 在线注释"模式，但免疫判定会降级。

## 示例数据说明

`examples/demo_sarscov2_timeseries.fasta` 是在 Wuhan-Hu-1 参考基因组上，按**真实密码子位置**工程化引入各关注株标志突变的
**合成演示序列**（非真实测序数据，序列头已标注 `SYNTHETIC DEMO`），用于离线跑通流程；突变清单见 `demo_mutation_log.tsv`。
实际研究请上传真实测序基因组。重新生成：`python examples/make_examples.py`。
`examples/report_*.html` 为预生成的示例分析报告。

## 已验证

- 变异调用对 30kb 基因组 ~10 ms；相邻/同聚区 SNP 不产生伪 indel；真实缺失/插入正确保留。
- 四时期演示样本的 S 蛋白突变与谱系签名（Delta/BA.1/XBB）判定正确。
- Streamlit 三页 `streamlit.testing` 无头测试全部通过（`test_app.py`）。
- UniProt（P0DTC2，110 个 feature）与 RCSB（2261 个 Spike 结构命中）在线接口实测可用。


## v0.2 新增：未来突变方向预测 + 三维结构影响

### ④ 未来突变方向预测（第④页签）
输出**可解释的重点监测突变排序（watchlist），不是确定性预言**。对 S 蛋白候选替换做六因子加权：

| 因子 | 含义 | 数据来源 |
|---|---|---|
| F1 免疫选择压力 | 是否中和/MHC表位、关键逃逸位点 | 知识库 |
| F2 结构可及性 | 残基相对溶剂可及面积 relSASA（表面才能被抗体接触） | 实验PDB 6M0J，Bio.PDB Shrake-Rupley |
| F3 受体界面 | 是否 RBD–ACE2 接触残基（重原子≤4.5Å） | 实验PDB 6M0J |
| F4 生化可容忍度 | BLOSUM62（结构是否可能允许） | 标准矩阵 |
| F5 流行率/逃逸先验 | 已在 Alpha/Delta/Omicron/XBB 观察到的突变、逃逸强度 | 知识快照 `sarscov2_prior.json`；有网可选 outbreak.info 实时频率，失败降级 |
| F6 时间轨迹 | 用户上传的多时期样本中晚期新出现/累积位点 | 第②页时间序列 |

输出 0–100 优先级、高/中/低置信度、方向分类（抗体/T细胞逃逸、受体结合/感染性、一般适应），
并逐条列出人类可读依据与六因子分值；支持内联 3Dmol.js 三维高亮（无网降级为表格）与 HTML 导出。

**已验证**：在 34 个候选中，14/14 个已知 RBD 逃逸位点（K417/K444/V445/G446/L452/F456/N460/T478/E484/F486/Q493/Q498/N501/Y505）全部被高优召回；
加入四时期轨迹后，晚期出现的 F486/K444 位点评分正确上调（+4.8/+4.8）。

### ③ 蛋白结构影响升级为真实三维评估
新增 `core/structure_impact.py`：下载并缓存实验 PDB（默认 6M0J RBD–ACE2），自动做 S 全长位点→PDB 残基映射，
逐突变计算 **relSASA（核心/部分埋藏/表面）、ACE2 界面接触、8Å 结构邻域、BLOSUM62**，给出结构后果叙述与"伪ΔΔG"失稳档。
实测与已知生物学一致：N501/Q498 为界面埋藏残基（relSASA<0.05），E484 为界面表面残基（0.48），D614 不在 RBD 覆盖区（正确标记为不可评估）。

> 方法边界：relSASA/界面为**实验结构实测**；"伪ΔΔG"是 BLOSUM×埋藏度启发式，**不是 FoldX/Rosetta 物理自由能**。
> 后续可接 FoldX/Rosetta ddG、ESMFold 突变体建模、netMHCpan 表位预测。实时频率依赖第三方 API 可用性。

## 当前局限与后续路线

1. 翻译使用标准遗传密码，暂未处理冠状病毒核糖体移码（ORF1ab）与 HIV 多重剪接/移码——重叠/移码区的氨基酸编号可能偏移；建议后续引入逐 ORF 的成熟注释。
2. 表位/逃逸判定依赖精选知识库覆盖；"结构严重性"是性质启发式，**不是 ΔΔG**。架构已预留重计算接口，下一步可接：
   - FoldX / Rosetta `ddg_monomer` 计算真实 ΔΔG；
   - ESMFold / AlphaFold 对突变体建模并做结构比对；
   - MHC 结合预测（netMHCpan）与表位全库（IEDB）覆盖更多病毒。
3. 谱系签名为简化集合匹配，生产环境建议对接 Pangolin / Nextclade。
4. 目前为单倍体一致性序列；尚未做测序 read 层面的低频变异（intrahost minor variant）分析。

## 目录

```
app.py                     Streamlit 入口
virinsight/core/          io/identify/variant/annotate/immuno/online/pipeline
virinsight/data/refs/      内置参考基因组
virinsight/data/kb/        病毒注册、基因与免疫表位知识库
virinsight/report.py       HTML 报告导出（含预测报告）
virinsight/data/pdb/       实验结构缓存(6M0J)
examples/                 合成演示序列、突变清单、示例报告
test_app.py               Streamlit 无头测试
requirements.txt
```
