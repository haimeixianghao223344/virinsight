"""在线数据库增强（可选，失败时优雅降级，绝不阻塞主流程）。
- UniProt REST：拉取活性位点 / 跨膜 / 结构域注释
- RCSB PDB：按 UniProt accession 查实验结构
返回结果统一标注来源与获取时间。
"""
from __future__ import annotations
import json
import datetime
import urllib.request
import urllib.error

UA = {"User-Agent": "VirInsight-prototype/0.1 (research; contact: local)"}
TIMEOUT = 12


def _get(url: str, accept: str = "application/json") -> str:
    req = urllib.request.Request(url, headers={**UA, "Accept": accept})
    with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
        return r.read().decode("utf-8")


def uniprot_features(accession: str) -> dict:
    """拉取 UniProt 关键 feature（Active site / Binding / Transmembrane / Domain / Region / Site）。"""
    out = {"ok": False, "accession": accession, "features": [], "error": None,
           "fetched_at": datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z",
           "source": "UniProt REST"}
    if not accession:
        out["error"] = "no accession"
        return out
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        data = json.loads(_get(url))
    except Exception as e:  # noqa: BLE001
        out["error"] = f"{type(e).__name__}: {e}"
        return out
    interesting = {"Active site", "Binding site", "Transmembrane", "Domain",
                   "Region", "Site", "Helix", "Disulfide bond", "Glycosylation site"}
    for f in data.get("features", []):
        if f.get("type") in interesting:
            loc = f.get("location", {})
            start = loc.get("start", {}).get("value")
            end = loc.get("end", {}).get("value")
            out["features"].append(dict(
                type=f.get("type"), start=start, end=end,
                description=f.get("description", ""),
            ))
    out["ok"] = True
    out["protein_name"] = (data.get("proteinDescription", {})
                           .get("recommendedName", {}).get("fullName", {}).get("value", ""))
    return out


def rcsb_by_uniprot(accession: str) -> dict:
    """查询 RCSB 中引用该 UniProt 的实验结构。"""
    out = {"ok": False, "accession": accession, "structures": [], "error": None,
           "fetched_at": datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z",
           "source": "RCSB PDB REST"}
    if not accession:
        out["error"] = "no accession"
        return out
    query = {
        "query": {"type": "terminal", "service": "text",
                  "parameters": {"attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                                 "operator": "exact_match", "value": accession}},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 15}}}
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    try:
        body = json.dumps(query).encode("utf-8")
        req = urllib.request.Request(url, data=body,
                                     headers={**UA, "Content-Type": "application/json",
                                              "Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
            data = json.loads(r.read().decode("utf-8"))
        for r in data.get("result_set", []):
            out["structures"].append(r["identifier"])
        out["ok"] = True
    except Exception as e:  # noqa: BLE001
        out["error"] = f"{type(e).__name__}: {e}"
    return out


def annotate_online(accession: str) -> dict:
    return {"uniprot": uniprot_features(accession), "rcsb": rcsb_by_uniprot(accession)}
import urllib.parse  # noqa: E402  (used in rcsb_by_uniprot)
