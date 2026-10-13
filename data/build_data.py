#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从 RCSB PDB 文件抽取 CA 主链 + 二级结构 + 血红素辅基，输出为 JS 数据文件。"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def parse_pdb(path, het_resnames=("HEM",)):
    chains = {}
    helices = []
    sheets = []
    hets = []
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec == "ATOM":
                name = line[12:16].strip()
                alt = line[16]
                if alt not in (" ", "A"):
                    continue
                if name != "CA":
                    continue
                chain = line[21]
                resname = line[17:20].strip()
                try:
                    resseq = int(line[22:26])
                except ValueError:
                    resseq = 0
                ch = chains.setdefault(chain, [])
                # 跳过同一残基重复
                if ch and ch[-1]["seq"] == resseq:
                    continue
                ch.append({
                    "seq": resseq,
                    "rn": resname,
                    "xyz": [round(float(line[30:38]), 3),
                            round(float(line[38:46]), 3),
                            round(float(line[46:54]), 3)],
                    "ss": "C",
                })
            elif rec == "HELIX":
                # initChainID 在第20列(idx19)，initSeqNum 第22-25列(idx21:25)，endSeqNum 第34-37列(idx33:37)
                try:
                    helices.append((line[19], int(line[21:25]), int(line[33:37])))
                except ValueError:
                    pass
            elif rec == "SHEET":
                try:
                    sheets.append((line[21], int(line[22:26]), int(line[33:37])))
                except ValueError:
                    pass
            elif rec == "HETATM":
                resname = line[17:20].strip()
                if resname not in het_resnames:
                    continue
                name = line[12:16].strip()
                elem = line[76:78].strip().upper()
                if not elem:
                    elem = name[0]
                if elem in ("H", "D"):
                    continue
                hets.append({
                    "chain": line[21],
                    "rn": resname,
                    "el": elem if elem != "FE" else "Fe",
                    "n": name,
                    "xyz": [round(float(line[30:38]), 3),
                            round(float(line[38:46]), 3),
                            round(float(line[46:54]), 3)],
                })
    # 标注二级结构
    for ch_id, residues in chains.items():
        for r in residues:
            s = r["seq"]
            for (c, a, b) in helices:
                if c == ch_id and a <= s <= b:
                    r["ss"] = "H"
            for (c, a, b) in sheets:
                if c == ch_id and a <= s <= b:
                    r["ss"] = "E"
    return chains, hets


def chain_info(chains, hets):
    out_chains = []
    for cid in sorted(chains):
        res = chains[cid]
        out_chains.append({
            "id": cid,
            "n": len(res),
            "first": res[0]["rn"] if res else "",
            "xyz": [r["xyz"] for r in res],
            "ss": "".join(r["ss"] for r in res),
        })
    out_hets = []
    for cid in sorted(chains):
        atoms = [h for h in hets if h["chain"] == cid]
        if not atoms:
            continue
        out_hets.append({
            "chain": cid,
            "atoms": [{"el": a["el"], "n": a["n"], "xyz": a["xyz"]} for a in atoms],
        })
    return out_chains, out_hets


def parse_segment(path, chain_id, lo, hi):
    """抽取一段连续肽链的主链原子（N, CA, C, O）与 CB，用于展示真实二级结构。"""
    out = {}
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[21] != chain_id:
                continue
            name = line[12:16].strip()
            if name not in ("N", "CA", "C", "O", "CB"):
                continue
            if line[16] not in (" ", "A"):
                continue
            try:
                seq = int(line[22:26])
            except ValueError:
                continue
            if not (lo <= seq <= hi):
                continue
            resname = line[17:20].strip()
            r = out.setdefault(seq, {"seq": seq, "rn": resname})
            r[name] = [round(float(line[30:38]), 3),
                       round(float(line[38:46]), 3),
                       round(float(line[46:54]), 3)]
    res = [out[k] for k in sorted(out)]
    # 去除缺原子的残基
    return [r for r in res if all(k in r for k in ("N", "CA", "C", "O"))]


def parse_ssbonds(path, keep_chains):
    """抽取指定链之间的二硫键（含 SG 硫原子坐标）。"""
    cys = {}
    bonds = []
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec == "ATOM":
                if line[21] not in keep_chains:
                    continue
                if line[12:16].strip() != "SG":
                    continue
                if line[16] not in (" ", "A"):
                    continue
                try:
                    seq = int(line[22:26])
                except ValueError:
                    continue
                cys[(line[21], seq)] = [round(float(line[30:38]), 3),
                                        round(float(line[38:46]), 3),
                                        round(float(line[46:54]), 3)]
            elif rec == "SSBOND":
                try:
                    c1, s1, c2, s2 = line[15], int(line[17:21]), line[29], int(line[31:35])
                except ValueError:
                    continue
                if c1 in keep_chains and c2 in keep_chains:
                    bonds.append({"a": [c1, s1], "b": [c2, s2]})
    out = []
    for b in bonds:
        pa, pb = cys.get(tuple(b["a"])), cys.get(tuple(b["b"]))
        if not pa or not pb:
            continue
        d = sum((pa[i] - pb[i]) ** 2 for i in range(3)) ** 0.5
        out.append({"a": b["a"][0] + str(b["a"][1]), "b": b["b"][0] + str(b["b"][1]),
                    "ca": pa, "cb": pb, "d": round(d, 2)})
    return cys, out


def main():
    data = {}
    # 真实二级结构片段：蛋白质G B1结构域（1GB1）
    segs = {}
    for key, (fname, cid, lo, hi, label, note) in {
        "helix_1gb1": ("1GB1.pdb", "A", 22, 36, "真实 α-螺旋", "1GB1 第22–36位残基"),
        "hairpin_1gb1": ("1GB1.pdb", "A", 42, 55, "真实反平行 β-发夹", "1GB1 第42–55位残基，两条β股反向平行"),
    }.items():
        res = parse_segment(os.path.join(HERE, fname), cid, lo, hi)
        segs[key] = {"label": label, "note": note, "pdb": "1GB1", "chain": cid, "residues": res}
        print(f"{key}: {len(res)} 残基 ({res[0]['seq']}-{res[-1]['seq']})")
    data["segments"] = segs
    for key, fname, label, desc in [
        ("myoglobin", "1MBO.pdb", "肌红蛋白", "抹香鲸肌红蛋白 · 1条肽链 · 153个氨基酸"),
        ("hemoglobin", "2HHB.pdb", "血红蛋白", "人血红蛋白 · 4条肽链（2α+2β）· 574个氨基酸"),
    ]:
        path = os.path.join(HERE, fname)
        chains, hets = parse_pdb(path)
        out_chains, out_hets = chain_info(chains, hets)
        data[key] = {
            "label": label,
            "pdb": key.upper() if key == "myoglobin" else "2HHB",
            "desc": desc,
            "chains": out_chains,
            "hets": out_hets,
        }
        n_aa = sum(c["n"] for c in out_chains)
        n_h = sum(len(h["atoms"]) for h in out_hets)
        print(f"{key}: {len(out_chains)} 链, {n_aa} 氨基酸, 辅基原子 {n_h}")

    # 胰岛素：A 链 21 个氨基酸 + B 链 30 个氨基酸，链间 2 个二硫键 + 链内 1 个
    cys, ssb = parse_ssbonds(os.path.join(HERE, "4INS.pdb"), "AB")
    chains, hets = parse_pdb(os.path.join(HERE, "4INS.pdb"))
    out_chains, _ = chain_info(chains, hets)
    out_chains = [c for c in out_chains if c["id"] in ("A", "B")]
    data["insulin"] = {
        "label": "胰岛素",
        "pdb": "4INS",
        "desc": "胰岛素 · 2条肽链（A链21个氨基酸 + B链30个氨基酸）· 2个链间二硫键",
        "chains": out_chains,
        "hets": [],
        "ssbonds": ssb,
        "cys": {k[0] + str(k[1]): v for k, v in cys.items()},
    }
    print("insulin:", [(c["id"], c["n"]) for c in out_chains], "二硫键", [(s["a"], s["b"], s["d"]) for s in ssb])

    js = "window.PROTEIN_PDB = " + json.dumps(data, ensure_ascii=False, separators=(",", ":")) + ";\n"
    with open(os.path.join(HERE, "proteins.js"), "w") as fh:
        fh.write(js)
    print("written proteins.js", round(len(js) / 1024, 1), "KB")


if __name__ == "__main__":
    main()
