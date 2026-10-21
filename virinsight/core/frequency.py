"""实时突变频率客户端（best-effort）。

环境内 outbreak.info / LAPIS(cov-spectrum) 可能不可达，因此：
- 任何失败都返回 ok=False，调用方使用本地先验，绝不阻塞；
- 仅对重点候选位点查询，控制请求数；
- 频率与增长方向作为先验之外的【实时证据】，单独标注来源与时间。
"""
from __future__ import annotations
import datetime
import json
import urllib.request
import urllib.parse

UA = {"User-Agent": "VirInsight/0.2 (research)", "Accept": "application/json"}
TIMEOUT = 12


def _get_json(url: str, headers=None):
    req = urllib.request.Request(url, headers={**UA, **(headers or {})})
    with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
        return json.loads(r.read().decode("utf-8"))


def outbreak_mutation_prevalence(mutation: str) -> dict:
    """outbreak.info 单点全球频率。mutation 形如 S:484K。"""
    out = {"ok": False, "source": "outbreak.info", "mutation": mutation,
           "fetched_at": datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z"}
    url = ("https://api.outbreak.info/genomics/mutation-details?mutations="
           + urllib.parse.quote(mutation))
    try:
        data = _get_json(url)
        if data and isinstance(data, dict):
            res = data.get(mutation) or next(iter(data.values()))
            if isinstance(res, dict):
                out["ok"] = True
                out["prevalence"] = res.get("prevalence")
                out["total_count"] = res.get("mutation_count")
    except Exception as e:  # noqa: BLE001
        out["error"] = f"{type(e).__name__}: {e}"
    return out


def fetch_live_prior(mutations: list[str], cap: int = 12) -> dict:
    """批量取实时频率，返回 {mutation: info}；整体失败返回空 dict。"""
    result = {}
    for m in mutations[:cap]:
        info = outbreak_mutation_prevalence(m)
        if info["ok"]:
            result[m] = info
    return result
