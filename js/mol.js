/* ============================================================
 * mol.js —— 分子几何构建与渲染（three.js r128）
 * 单位统一使用埃（Å），与 PDB 坐标一致
 * ============================================================ */
(function (global) {
  'use strict';
  const T = global.THREE;
  const V3 = T.Vector3;
  const D2R = Math.PI / 180;

  /* ---------------- 颜色与半径 ---------------- */
  const CPK = {
    C: 0x515a66, N: 0x2c6fd4, O: 0xdc4b3a, H: 0xdde3ea,
    S: 0xd9a406, Fe: 0xc2452b, P: 0xe07a2f,
  };
  const RAD = { C: 0.40, N: 0.38, O: 0.36, H: 0.245, S: 0.46, Fe: 0.46, P: 0.42 };
  const SS_COLOR = { H: 0xe8703f, E: 0x3f8fd6, C: 0xa9b4c0 };

  function lighten(hex, amt) {
    const c = new T.Color(hex);
    c.lerp(new T.Color(0xffffff), amt);
    return c.getHex();
  }

  /* ---------------- 几何工具 ---------------- */
  // 由 A-B-C 三个已知原子，按 (键长 r, 键角 BCD=theta, 二面角 ABCD=phi) 放置 D
  function placeAtom(A, B, C, r, thetaDeg, phiDeg) {
    const bcDir = C.clone().sub(B).normalize();
    const abDir = B.clone().sub(A);
    let n = new V3().crossVectors(bcDir, abDir);
    if (n.lengthSq() < 1e-8) n.set(0.35, 0.7, 0.1).cross(bcDir);
    n.normalize();
    const m = new V3().crossVectors(n, bcDir).normalize();
    const t = thetaDeg * D2R, p = phiDeg * D2R;
    const e1 = bcDir.clone().negate(); // C -> B
    return new V3()
      .addScaledVector(e1, r * Math.cos(t))
      .addScaledVector(m, -r * Math.cos(p) * Math.sin(t))
      .addScaledVector(n, r * Math.sin(p) * Math.sin(t))
      .add(C);
  }

  function dihedral(a, b, c, d) {
    const b1 = b.clone().sub(a), b2 = c.clone().sub(b), b3 = d.clone().sub(c);
    const n1 = new V3().crossVectors(b1, b2).normalize();
    const n2 = new V3().crossVectors(b2, b3).normalize();
    const m1 = new V3().crossVectors(n1, b2.clone().normalize());
    return Math.atan2(m1.dot(n2), n1.dot(n2)) / D2R;
  }

  /* ---------------- 侧链（R 基）定义 ----------------
   * 第一个原子一定是 CB，用角平分线法放置（只依赖 N、CA、C）
   * 后续原子用内坐标：ref = [A, B, C] 三个参考原子下标（-1=CA, -2=N, 0=CB, 1=CG...）
   */
  const SIDECHAINS = {
    Gly: { cn: '甘氨酸', abbr: 'Gly', r: '—H', type: '非极性', atoms: [], bonds: [] },
    Ala: { cn: '丙氨酸', abbr: 'Ala', r: '—CH₃', type: '非极性', atoms: [{ el: 'C', n: 'CB' }], bonds: [] },
    Val: {
      cn: '缬氨酸', abbr: 'Val', r: '—CH(CH₃)₂', type: '非极性（必需）',
      atoms: [
        { el: 'C', n: 'CB' },
        { el: 'C', n: 'CG1', ref: [-2, -1, 0], r: 1.53, a: 111, d: 60 },
        { el: 'C', n: 'CG2', ref: [-2, -1, 0], r: 1.53, a: 111, d: -60 },
      ],
      bonds: [[0, 1], [0, 2]],
    },
    Leu: {
      cn: '亮氨酸', abbr: 'Leu', r: '—CH₂CH(CH₃)₂', type: '非极性（必需）',
      atoms: [
        { el: 'C', n: 'CB' },
        { el: 'C', n: 'CG', ref: [-2, -1, 0], r: 1.53, a: 112, d: 180 },
        { el: 'C', n: 'CD1', ref: [-1, 0, 1], r: 1.53, a: 112, d: 60 },
        { el: 'C', n: 'CD2', ref: [-1, 0, 1], r: 1.53, a: 112, d: -60 },
      ],
      bonds: [[0, 1], [1, 2], [1, 3]],
    },
    Ser: {
      cn: '丝氨酸', abbr: 'Ser', r: '—CH₂OH', type: '极性不带电',
      atoms: [
        { el: 'C', n: 'CB' },
        { el: 'O', n: 'OG', ref: [-2, -1, 0], r: 1.42, a: 110.8, d: 180 },
      ],
      bonds: [[0, 1]],
    },
    Thr: {
      cn: '苏氨酸', abbr: 'Thr', r: '—CH(OH)CH₃', type: '极性不带电（必需）',
      atoms: [
        { el: 'C', n: 'CB' },
        { el: 'O', n: 'OG1', ref: [-2, -1, 0], r: 1.43, a: 110, d: 60 },
        { el: 'C', n: 'CG2', ref: [-2, -1, 0], r: 1.53, a: 111, d: -60 },
      ],
      bonds: [[0, 1], [0, 2]],
    },
    Cys: {
      cn: '半胱氨酸', abbr: 'Cys', r: '—CH₂SH', type: '极性（含硫）',
      atoms: [
        { el: 'C', n: 'CB' },
        { el: 'S', n: 'SG', ref: [-2, -1, 0], r: 1.81, a: 113.8, d: 180 },
      ],
      bonds: [[0, 1]],
    },
    Lys: {
      cn: '赖氨酸', abbr: 'Lys', r: '—(CH₂)₄NH₂', type: '碱性（必需）',
      atoms: [
        { el: 'C', n: 'CB' },
        { el: 'C', n: 'CG', ref: [-2, -1, 0], r: 1.53, a: 112, d: 180 },
        { el: 'C', n: 'CD', ref: [-1, 0, 1], r: 1.53, a: 112, d: 180 },
        { el: 'C', n: 'CE', ref: [0, 1, 2], r: 1.53, a: 112, d: 180 },
        { el: 'N', n: 'NZ', ref: [1, 2, 3], r: 1.49, a: 111, d: 180 },
      ],
      bonds: [[0, 1], [1, 2], [2, 3], [3, 4]],
    },
    Glu: {
      cn: '谷氨酸', abbr: 'Glu', r: '—CH₂CH₂COOH', type: '酸性',
      atoms: [
        { el: 'C', n: 'CB' },
        { el: 'C', n: 'CG', ref: [-2, -1, 0], r: 1.53, a: 112, d: 180 },
        { el: 'C', n: 'CD', ref: [-1, 0, 1], r: 1.53, a: 112, d: 180 },
        { el: 'O', n: 'OE1', ref: [0, 1, 2], r: 1.25, a: 118, d: 0 },
        { el: 'O', n: 'OE2', ref: [0, 1, 2], r: 1.25, a: 118, d: 180 },
      ],
      bonds: [[0, 1], [1, 2], [2, 3], [2, 4]],
    },
    Phe: {
      cn: '苯丙氨酸', abbr: 'Phe', r: '—CH₂—（苯环）', type: '非极性（芳香族，必需）',
      atoms: [
        { el: 'C', n: 'CB' },
        { el: 'C', n: 'CG', ref: [-2, -1, 0], r: 1.51, a: 113, d: 90 },
        { el: 'C', n: 'CD1', ref: [-1, 0, 1], r: 1.39, a: 120, d: 0 },
        { el: 'C', n: 'CD2', ref: [-1, 0, 1], r: 1.39, a: 120, d: 180 },
        { el: 'C', n: 'CE1', ref: [0, 1, 2], r: 1.39, a: 120, d: 180 },
        { el: 'C', n: 'CE2', ref: [0, 1, 3], r: 1.39, a: 120, d: 180 },
        { el: 'C', n: 'CZ', ref: [1, 2, 4], r: 1.39, a: 120, d: 0 },
      ],
      bonds: [[0, 1], [1, 2], [1, 3], [2, 4], [3, 5], [4, 6], [5, 6]],
    },
  };

  const DEFAULT_R = ['Gly', 'Ala', 'Val', 'Ser', 'Cys', 'Lys', 'Glu', 'Phe'];

  /* ---------------- 肽链数据构建 ---------------- */
  function buildPeptideAtoms(opts) {
    opts = opts || {};
    const n = Math.max(1, opts.n || 1);
    const phi = opts.phi === undefined ? -120 : opts.phi;
    const psi = opts.psi === undefined ? 120 : opts.psi;
    const scList = opts.sidechains || [];
    const showH = opts.showH !== false;
    const showO = opts.showO !== false;
    const termini = opts.termini !== false;

    const atoms = [];
    const bonds = [];
    const residues = [];
    const push = (el, name, res, pos, tag) => {
      atoms.push({ el, name, res, pos, tag });
      return atoms.length - 1;
    };
    const bond = (i, j, order, kind) => {
      if (i === undefined || j === undefined || i < 0 || j < 0) return;
      bonds.push({ a: i, b: j, order: order || 1, kind: kind || 'single' });
    };

    // ---- 种子：N1, CA1, C1（键角 N-CA-C = 111°）----
    let N = new V3(0, 0, 0);
    let CA = new V3(1.458, 0, 0);
    let C = new V3(1.458 + 1.525 * Math.cos(69 * D2R), 1.525 * Math.sin(69 * D2R), 0);
    const Ns = [N], CAs = [CA], Cs = [C];
    for (let i = 1; i < n; i++) {
      const Nn = placeAtom(Ns[i - 1], CAs[i - 1], Cs[i - 1], 1.329, 116, psi);
      const CAn = placeAtom(CAs[i - 1], Cs[i - 1], Nn, 1.458, 121, 180);
      const Cn = placeAtom(Cs[i - 1], Nn, CAn, 1.525, 111, phi);
      Ns.push(Nn); CAs.push(CAn); Cs.push(Cn);
    }

    // ---- 羰基 O ----
    const Os = [];
    for (let i = 0; i < n; i++) {
      Os.push(placeAtom(Ns[i], CAs[i], Cs[i], 1.231, 121, psi + 180));
    }

    // ---- 逐个残基登记原子 ----
    for (let i = 0; i < n; i++) {
      const rec = { idx: i };
      rec.N = push('N', 'N', i, Ns[i], 'bb');
      rec.CA = push('C', 'CA', i, CAs[i], 'bb');
      rec.C = push('C', 'C', i, Cs[i], 'bb');
      if (showO) rec.O = push('O', 'O', i, Os[i], 'bb');

      // 主链键
      bond(rec.N, rec.CA);
      bond(rec.CA, rec.C);
      if (rec.O !== undefined) bond(rec.C, rec.O, 2, 'co');
      if (i > 0) {
        const prev = residues[i - 1];
        bond(prev.C, rec.N, 1, 'peptide');
      }

      // 氨基上的 H
      rec.HN = [];
      if (showH) {
        if (i === 0) {
          const dummy = Ns[0].clone().addScaledVector(Ns[0].clone().sub(CAs[0]).normalize(), 1.0);
          rec.HN.push(push('H', 'H1', i, placeAtom(dummy, CAs[i], Ns[i], 1.01, 110, 120), 'hn'));
          rec.HN.push(push('H', 'H2', i, placeAtom(dummy, CAs[i], Ns[i], 1.01, 110, -120), 'hn'));
          rec.HN.forEach((h) => bond(rec.N, h, 1, 'hn'));
        } else {
          // 反式肽键：酰胺 H 与前一个羰基 O 处于肽平面内、分居 C—N 键两侧（二面角 180°）
          const h = push('H', 'H', i, placeAtom(Os[i - 1], Cs[i - 1], Ns[i], 1.01, 119, 180), 'hn');
          rec.HN.push(h);
          bond(rec.N, h, 1, 'hn');
        }
      }

      // α-氢
      if (showH && opts.hAlpha !== false) {
        const nv = Ns[i].clone().sub(CAs[i]).normalize();
        const cv = Cs[i].clone().sub(CAs[i]).normalize();
        const axis = nv.clone().add(cv).normalize();
        const perp = new V3().crossVectors(nv, cv).normalize();
        const h = CAs[i].clone()
          .addScaledVector(axis, -1.09 * Math.cos(55.5 * D2R))
          .addScaledVector(perp, -1.09 * Math.sin(55.5 * D2R));
        rec.HA = push('H', 'HA', i, h, 'ha');
        bond(rec.CA, rec.HA, 1, 'ha');
      }

      // 羧基末端 —OH
      if (termini && i === n - 1 && showH) {
        const dummy = Cs[i].clone().addScaledVector(Cs[i].clone().sub(CAs[i]).normalize(), 1.0);
        const dO = dihedral(dummy, CAs[i], Cs[i], Os[i]);
        // 羧基与 CA 共面，O 与 OXT 在该平面内分居 C—CA 键两侧（二面角相差 180°）
        rec.OXT = push('O', 'OXT', i, placeAtom(dummy, CAs[i], Cs[i], 1.34, 116, dO + 180), 'cooh');
        rec.HXT = push('H', 'HXT', i, placeAtom(Os[i], Cs[i], rec.OXT ? atoms[rec.OXT].pos : Os[i], 0.98, 107, 0), 'cooh');
        bond(rec.C, rec.OXT, 1, 'cooh');
        bond(rec.OXT, rec.HXT, 1, 'cooh');
      }
      residues.push(rec);
    }

    // ---- 侧链（R 基）----
    if (opts.sidechains !== false) {
      for (let i = 0; i < n; i++) {
        const key = scList[i] || scList[0] || 'Ala';
        const def = SIDECHAINS[key] || SIDECHAINS.Ala;
        const rec = residues[i];
        const nv = Ns[i].clone().sub(CAs[i]).normalize();
        const cv = Cs[i].clone().sub(CAs[i]).normalize();
        const axis = nv.clone().add(cv).normalize();
        const perp = new V3().crossVectors(nv, cv).normalize();
        const cb = CAs[i].clone()
          .addScaledVector(axis, -1.53 * Math.cos(55.5 * D2R))
          .addScaledVector(perp, 1.53 * Math.sin(55.5 * D2R));
        const built = [];
        built.push(push('C', 'CB', i, cb, 'sc'));
        bond(rec.CA, built[0], 1, 'sc');
        (def.atoms || []).forEach((a, k) => {
          if (k === 0) return; // CB 已放置
          const ref = a.ref.map((x) => (x === -1 ? CAs[i] : x === -2 ? Ns[i] : atoms[built[x]].pos));
          built.push(push(a.el, a.n, i, placeAtom(ref[0], ref[1], ref[2], a.r, a.a, a.d), 'sc'));
        });
        (def.bonds || []).forEach(([x, y]) => bond(built[x], built[y], 1, 'sc'));
        rec.CB = built[0];
        rec.sc = built;
        rec.scKey = key;
      }
    }

    return { atoms, bonds, residues, Ns, CAs, Cs, Os, n };
  }

  /* ---------------- 渲染：球棍模型 ---------------- */
  const GEO_SPHERE = new T.SphereGeometry(1, 20, 14);
  const GEO_CYL = new T.CylinderGeometry(1, 1, 1, 12, 1, true);
  const YAXIS = new V3(0, 1, 0);
  const matCache = new Map();
  function mat(color, opts) {
    opts = opts || {};
    const key = color + '|' + (opts.shininess || 60) + '|' + (opts.opacity || 1);
    if (!opts.transparent && matCache.has(key)) return matCache.get(key);
    const m = new T.MeshPhongMaterial({
      color,
      shininess: opts.shininess === undefined ? 60 : opts.shininess,
      specular: 0x2a3340,
      transparent: !!opts.transparent,
      opacity: opts.opacity === undefined ? 1 : opts.opacity,
      side: opts.side || T.FrontSide,
    });
    if (!opts.transparent) matCache.set(key, m);
    return m;
  }

  function addSphere(parent, el, pos, radius, opts) {
    opts = opts || {};
    const m = new T.Mesh(GEO_SPHERE, opts.material || mat(opts.color === undefined ? CPK[el] : opts.color, opts));
    m.scale.setScalar(radius === undefined ? RAD[el] : radius);
    m.position.copy(pos);
    parent.add(m);
    return m;
  }

  function cylBetween(parent, pa, pb, radius, color, opts) {
    const dir = pb.clone().sub(pa);
    const len = dir.length();
    if (len < 1e-6) return null;
    const m = new T.Mesh(GEO_CYL, opts && opts.material ? opts.material : mat(color, opts));
    m.scale.set(radius, len, radius);
    m.position.copy(pa).addScaledVector(dir, 0.5);
    m.quaternion.setFromUnitVectors(YAXIS, dir.clone().normalize());
    parent.add(m);
    return m;
  }

  function addBond(parent, pa, pb, opts) {
    opts = opts || {};
    const radius = opts.radius || 0.115;
    const dir = pb.clone().sub(pa);
    if (dir.length() < 1e-6) return null;
    if (opts.order === 2) {
      let perp = new V3(0.3, 0.8, 0.15).cross(dir);
      if (perp.lengthSq() < 1e-6) perp = new V3(1, 0, 0).cross(dir);
      perp.normalize().multiplyScalar(0.105);
      const g = new T.Group();
      cylBetween(g, pa.clone().add(perp), pb.clone().add(perp), radius * 0.62, opts.colorB || 0xc3ccd6);
      cylBetween(g, pa.clone().sub(perp), pb.clone().sub(perp), radius * 0.62, opts.colorB || 0xc3ccd6);
      parent.add(g);
      return g;
    }
    const mid = pa.clone().addScaledVector(dir, 0.5);
    const g = new T.Group();
    cylBetween(g, pa, mid, radius, opts.colorA === undefined ? 0xb9c3ce : opts.colorA);
    cylBetween(g, mid, pb, radius, opts.colorB === undefined ? 0xb9c3ce : opts.colorB);
    parent.add(g);
    return g;
  }

  // 动态键：两端跟随两个 mesh 的位置（自动换算到 parent 的局部坐标系）
  function dynamicBond(parent, meshA, meshB, radius, color) {
    const m = new T.Mesh(GEO_CYL, mat(color || 0xb9c3ce));
    parent.add(m);
    const tmpA = new V3(), tmpB = new V3();
    const api = {
      mesh: m,
      radius: radius,
      visible: true,
      update() {
        if (!api.visible) { m.visible = false; return; }
        meshA.getWorldPosition(tmpA);
        meshB.getWorldPosition(tmpB);
        parent.worldToLocal(tmpA);
        parent.worldToLocal(tmpB);
        const dir = tmpB.clone().sub(tmpA);
        const len = dir.length();
        if (len < 1e-5) { m.visible = false; return; }
        m.visible = true;
        m.scale.set(api.radius, len, api.radius);
        m.position.copy(tmpA).addScaledVector(dir, 0.5);
        m.quaternion.setFromUnitVectors(YAXIS, dir.normalize());
      },
    };
    return api;
  }

  function dashedLine(parent, pa, pb, color, width) {
    const g = new T.BufferGeometry().setFromPoints([pa.clone(), pb.clone()]);
    const m = new T.LineDashedMaterial({ color: color || 0x1f8a70, dashSize: 0.34, gapSize: 0.26, linewidth: width || 1 });
    const line = new T.Line(g, m);
    line.computeLineDistances();
    parent.add(line);
    return line;
  }

  /* ---------------- 完整分子渲染 ---------------- */
  function renderMol(data, opts) {
    opts = opts || {};
    const root = new T.Group();
    const atomMeshes = data.atoms.map(() => null);
    const bondMeshes = [];
    const resGroups = [];
    for (let i = 0; i < data.n; i++) {
      const g = new T.Group();
      resGroups.push(g);
      root.add(g);
    }

    data.atoms.forEach((a, i) => {
      let parent = resGroups[a.res];
      if (a.detached || (opts.detached && opts.detached.indexOf(i) >= 0)) parent = root;
      const m = addSphere(parent, a.el, a.pos, opts.scale ? RAD[a.el] * opts.scale : RAD[a.el], { color: a.color });
      atomMeshes[i] = m;
      a.mesh = m;
    });

    data.bonds.forEach((b) => {
      const aa = data.atoms[b.a], bb = data.atoms[b.b];
      const pa = aa.pos, pb = bb.pos;
      const sameRes = aa.res === bb.res;
      let parent = sameRes ? resGroups[aa.res] : root;
      if (aa.detached || bb.detached) parent = root;
      const colorA = lighten(aa.color === undefined ? CPK[aa.el] : aa.color, 0.35);
      const colorB = lighten(bb.color === undefined ? CPK[bb.el] : bb.color, 0.35);
      const mesh = addBond(parent, pa, pb, {
        order: b.order, kind: b.kind, colorA, colorB,
        radius: (b.order === 2 ? 0.11 : 0.115) * (opts.scale || 1),
      });
      b.mesh = mesh;
      bondMeshes.push({ bond: b, mesh });
    });

    return { root, atomMeshes, bondMeshes, resGroups, data };
  }

  /* ---------------- 自动成键（用于辅基等无连接信息的原子团） ---------------- */
  function autoBonds(atoms) {
    const out = [];
    for (let i = 0; i < atoms.length; i++) {
      for (let j = i + 1; j < atoms.length; j++) {
        const ei = atoms[i].el, ej = atoms[j].el;
        const lim = (ei === 'Fe' || ej === 'Fe') ? 2.25 : 1.78;
        if (atoms[i].pos.distanceTo(atoms[j].pos) < lim) out.push([i, j]);
      }
    }
    return out;
  }

  function renderHet(parent, atomList, opts) {
    opts = opts || {};
    const pts = atomList.map((a) => ({ el: a.el, pos: new V3(a.xyz[0], a.xyz[1], a.xyz[2]) }));
    const meshes = pts.map((p) => addSphere(parent, p.el, p.pos, RAD[p.el] * (opts.scale || 1), { color: p.el === 'Fe' ? 0xc2452b : undefined }));
    autoBonds(pts).forEach(([i, j]) => {
      addBond(parent, pts[i].pos, pts[j].pos, {
        colorA: lighten(CPK[pts[i].el], 0.35),
        colorB: lighten(CPK[pts[j].el], 0.35),
        radius: 0.1 * (opts.scale || 1),
      });
    });
    return meshes;
  }

  /* ---------------- 蛋白质主链卡通（管状） ---------------- */
  // TubeGeometry 的顶点数 = (管段数+1) × (横截面段数+1)，颜色需按横截面展开
  function ssColors(ssStr, segs, radial) {
    radial = radial || 8;
    const n = ssStr.length;
    const per = radial + 1;
    const arr = new Float32Array((segs + 1) * per * 3);
    const c = new T.Color();
    for (let i = 0; i <= segs; i++) {
      const f = (i / segs) * (n - 1);
      const idx = Math.min(n - 1, Math.max(0, Math.round(f)));
      c.set(SS_COLOR[ssStr[idx]] || SS_COLOR.C);
      for (let j = 0; j < per; j++) {
        const k = (i * per + j) * 3;
        arr[k] = c.r; arr[k + 1] = c.g; arr[k + 2] = c.b;
      }
    }
    return arr;
  }

  class Cartoon {
    constructor(points, opts) {
      opts = opts || {};
      this.points = points.map((p) => p.clone());
      this.radius = opts.radius || 0.85;
      this.segPerRes = opts.segPerRes || 6;
      this.group = new T.Group();
      this.color = opts.color;
      this.opacity = opts.opacity === undefined ? 1 : opts.opacity;
      const segs = Math.max(8, (this.points.length - 1) * this.segPerRes);
      this.segs = segs;
      if (this.color === undefined) {
        this.colors = ssColors(opts.ss || 'C'.repeat(this.points.length), segs);
      }
      this.material = new T.MeshPhongMaterial({
        vertexColors: this.color === undefined,
        color: this.color === undefined ? 0xffffff : this.color,
        shininess: 55,
        specular: 0x33404d,
        transparent: this.opacity < 1,
        opacity: this.opacity,
      });
      this.mesh = null;
      this.rebuild();
    }
    _geom() {
      const curve = new T.CatmullRomCurve3(this.points, false, 'centripetal');
      return new T.TubeGeometry(curve, this.segs, this.radius, 8, false);
    }
    rebuild() {
      const g = this._geom();
      if (this.colors) g.setAttribute('color', new T.BufferAttribute(this.colors, 3));
      if (this.mesh) { this.group.remove(this.mesh); this.mesh.geometry.dispose(); }
      this.mesh = new T.Mesh(g, this.material);
      this.group.add(this.mesh);
    }
    setPoints(pts) {
      for (let i = 0; i < this.points.length; i++) this.points[i].copy(pts[i]);
      this.rebuild();
    }
    dispose() {
      if (this.mesh) this.mesh.geometry.dispose();
      this.material.dispose();
    }
  }

  /* ---------------- 刚体叠合：把 src 三点组贴合到 dst 三点组 ---------------- */
  function basisFrom3(p0, p1, p2) {
    const e1 = p1.clone().sub(p0).normalize();
    const e3 = new V3().crossVectors(e1, p2.clone().sub(p0));
    if (e3.lengthSq() < 1e-8) e3.set(0, 0, 1).cross(e1);
    e3.normalize();
    const e2 = new V3().crossVectors(e3, e1).normalize();
    return new T.Matrix4().makeBasis(e1, e2, e3);
  }

  // 返回 {q, t}：p' = q * p + t
  function rigidFit(src, dst) {
    const cs = centerPoints(src).clone();
    const cd = centerPoints(dst).clone();
    const Ms = basisFrom3(src[0], src[1], src[2]);
    const Md = basisFrom3(dst[0], dst[1], dst[2]);
    const R = Md.multiply(Ms.clone().transpose());
    const q = new T.Quaternion().setFromRotationMatrix(R);
    const t = cd.clone().sub(cs.clone().applyQuaternion(q));
    return { q, t };
  }

  function applyRigid(pts, fit) {
    return pts.map((p) => p.clone().applyQuaternion(fit.q).add(fit.t));
  }

  function centerPoints(pts) {
    const c = new V3();
    pts.forEach((p) => c.add(p));
    c.divideScalar(pts.length || 1);
    return c;
  }

  /* ---------------- 随机卷曲（用于变性演示） ---------------- */
  function mulberry32(seed) {
    return function () {
      seed |= 0; seed = (seed + 0x6d2b79f5) | 0;
      let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  // 生成一段"去折叠"的主链 CA 轨迹（键长键角严格正确，二面角随机）
  function randomCoilCA(n, seed) {
    const rnd = mulberry32(seed || 20261013);
    const data = [];
    let phi = -120, psi = 120;
    // 用逐残基随机的 phi/psi 生成主链，再取 CA
    let N = new V3(0, 0, 0);
    let CA = new V3(1.458, 0, 0);
    let C = new V3(1.458 + 1.525 * Math.cos(69 * D2R), 1.525 * Math.sin(69 * D2R), 0);
    data.push(CA.clone());
    for (let i = 1; i < n; i++) {
      phi = -150 + rnd() * 100;      // -150 ~ -50
      psi = 95 + rnd() * 90;         // 95 ~ 185
      if (rnd() < 0.12) { psi = -40 + rnd() * 40; phi = -70 + rnd() * 25; }
      const Nn = placeAtom(N, CA, C, 1.329, 116, psi);
      const CAn = placeAtom(CA, C, Nn, 1.458, 121, 180);
      const Cn = placeAtom(C, Nn, CAn, 1.525, 111, phi);
      N = Nn; CA = CAn; C = Cn;
      data.push(CA.clone());
    }
    return data;
  }

  /* ---------------- 清理 ---------------- */
  function disposeGroup(group) {
    group.traverse((o) => {
      if (o.geometry && o.geometry !== GEO_SPHERE && o.geometry !== GEO_CYL) o.geometry.dispose();
      if (o.material && o.material.__shared !== true) {
        // 共享材质在缓存中，不销毁
      }
    });
  }

  global.MOL = {
    CPK, RAD, SS_COLOR, SIDECHAINS, DEFAULT_R,
    placeAtom, dihedral, lighten, mat,
    buildPeptideAtoms, renderMol, renderHet, autoBonds,
    addSphere, addBond, cylBetween, dynamicBond, dashedLine,
    Cartoon, centerPoints, randomCoilCA, mulberry32,
    rigidFit, applyRigid, basisFrom3,
  };
})(window);
