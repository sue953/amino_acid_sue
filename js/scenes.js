/* ============================================================
 * scenes.js —— 8 个教学场景
 * 沪科版2020 必修1 第2章第2节《蛋白质和核酸是重要的生物大分子》
 * ============================================================ */
(function (global) {
  'use strict';
  const T = global.THREE;
  const V3 = T.Vector3;
  const MOL = global.MOL;

  /* ---------------- 通用工具 ---------------- */
  function centerData(data) {
    const c = MOL.centerPoints(data.atoms.map((a) => a.pos));
    data.atoms.forEach((a) => a.pos.sub(c));
    return c;
  }
  function smoothstep(a, b, x) {
    if (b === a) return x < a ? 0 : 1;
    const t = Math.min(1, Math.max(0, (x - a) / (b - a)));
    return t * t * (3 - 2 * t);
  }
  function outward(p, center, d) {
    const v = p.clone().sub(center);
    if (v.lengthSq() < 1e-6) v.set(0, 1, 0);
    return p.clone().addScaledVector(v.normalize(), d);
  }
  // 由 PDB 片段构建可渲染的原子/键数据
  function buildSegment(residues, opts) {
    opts = opts || {};
    const atoms = [], bonds = [], res = [];
    const v = (a) => new V3(a[0], a[1], a[2]);
    const push = (el, name, i, pos, tag) => { atoms.push({ el, name, res: i, pos, tag }); return atoms.length - 1; };
    const bd = (i, j, order, kind) => bonds.push({ a: i, b: j, order: order || 1, kind: kind || 'single' });
    residues.forEach((r, i) => {
      const rec = { idx: i, rn: r.rn, seq: r.seq };
      rec.N = push('N', 'N', i, v(r.N), 'bb');
      rec.CA = push('C', 'CA', i, v(r.CA), 'bb');
      rec.C = push('C', 'C', i, v(r.C), 'bb');
      rec.O = push('O', 'O', i, v(r.O), 'bb');
      bd(rec.N, rec.CA); bd(rec.CA, rec.C); bd(rec.C, rec.O, 2, 'co');
      if (i > 0) bd(res[i - 1].C, rec.N, 1, 'peptide');
      if (opts.hAmide !== false && i > 0) {
        const p = residues[i - 1];
        rec.H = push('H', 'H', i,
          MOL.placeAtom(v(p.O), v(p.C), v(r.N), 1.01, 119, 180), 'hn');
        bd(rec.N, rec.H, 1, 'hn');
      }
      if (r.CB && opts.showCB !== false) {
        rec.CB = push('C', 'CB', i, v(r.CB), 'sc');
        bd(rec.CA, rec.CB, 1, 'sc');
      }
      res.push(rec);
    });
    return { atoms, bonds, residues: res, n: res.length };
  }
  // 主链氢键自动识别：H(i) 与 O(j)，|i-j| >= sep，距离 < 2.5 Å
  function findHBonds(data, sep) {
    const H = [], O = [], out = [];
    data.residues.forEach((r, i) => {
      if (r.H !== undefined) H.push({ res: i, idx: r.H });
      O.push({ res: i, idx: r.O });
    });
    H.forEach((h) => O.forEach((o) => {
      if (Math.abs(h.res - o.res) < sep) return;
      const d = data.atoms[h.idx].pos.distanceTo(data.atoms[o.idx].pos);
      if (d < 2.5) out.push({ h: h.idx, o: o.idx, d });
    }));
    return out;
  }
  function drawHBonds(parent, data, hb, color) {
    hb.forEach((b) => {
      MOL.dashedLine(parent, data.atoms[b.h].pos, data.atoms[b.o].pos, color || 0x1f8a70);
    });
  }
  function pdbChain(key, ci) {
    return global.PROTEIN_PDB[key].chains[ci];
  }
  function chainPoints(ch) {
    return ch.xyz.map((p) => new V3(p[0], p[1], p[2]));
  }

  /* ============================================================
   * ① 氨基酸
   * ============================================================ */
  const s1 = {
    id: 'amino',
    nav: '① 氨基酸',
    title: '氨基酸 —— 蛋白质的基本单位',
    lead: '蛋白质由氨基酸组成。组成生物体蛋白质的氨基酸已知有 22 种，其中常见的有 20 种。',
    points: [
      '每种氨基酸分子<b>至少</b>含有一个氨基（—NH₂）和一个羧基（—COOH）。',
      '并且都有<b>一个氨基和一个羧基连接在同一个碳原子上</b>，这个碳原子叫 α-碳原子（中心碳原子）。',
      'α-碳原子上还连接着<b>一个氢原子</b>和一个<b>侧链 R 基</b>。',
      '结构通式：<b>H₂N—CH(R)—COOH</b>。',
    ],
    exam: '判断某化合物是不是构成蛋白质的氨基酸，关键看：<b>是否有一个 —NH₂ 和一个 —COOH 连在同一个碳原子上</b>。R 基中额外出现的氨基或羧基不算。',
    build(ctx) {
      const g = new T.Group(); ctx.root.add(g);
      const data = MOL.buildPeptideAtoms({ n: 1, sidechains: ['Ala'], termini: true });
      centerData(data);
      const mol = MOL.renderMol(data, { scale: 1.4 });
      g.add(mol.root);
      const R = data.residues[0];
      const c0 = new V3(0, 0, 0);

      const L = [];
      L.push(ctx.lab('氨基 —NH₂', outward(data.atoms[R.N].pos, c0, 2.2), 'c-blue'));
      L.push(ctx.lab('羧基 —COOH', outward(data.atoms[R.C].pos, c0, 2.4), 'c-red'));
      L.push(ctx.lab('α-碳原子<br><span class="s">（中心碳原子）</span>', outward(data.atoms[R.CA].pos, c0, 2.0), 'c-ink'));
      L.push(ctx.lab('氢原子', outward(data.atoms[R.HA].pos, c0, 1.6), 'c-gray'));
      L.push(ctx.lab('R 基（侧链）', outward(data.atoms[R.CB].pos, c0, 2.0), 'c-teal'));
      ctx.setView(new V3(0, 0, 0), 12.5, new V3(0.15, 0.25, 1));
      ctx.stat('<div class="formula">H₂N — CH(R) — COOH</div><div class="note">拖动模型可以从不同角度观察四个基团的连接关系</div>');
      return {
        dispose() { L.forEach((l) => l.remove()); },
      };
    },
  };

  /* ============================================================
   * ② R 基决定种类
   * ============================================================ */
  const AA_TYPES = [
    ['Gly', '#8aa0b4'], ['Ala', '#5eb0a8'], ['Val', '#7b8ce0'], ['Ser', '#e0a63f'],
    ['Cys', '#d9d15a'], ['Lys', '#6a8de0'], ['Glu', '#e07a5f'], ['Phe', '#b07ad9'],
  ];
  const s2 = {
    id: 'rgroup',
    nav: '② R 基与种类',
    title: 'R 基不同 —— 氨基酸的种类不同',
    lead: '20 种常见氨基酸的区别，只在于 R 基（侧链）不同。R 基的化学性质决定了氨基酸的极性、酸碱性。',
    points: [
      '不同氨基酸的<b>主链部分完全相同</b>，唯一区别是 <b>R 基</b>不同，因此 R 基决定了氨基酸的种类。',
      '按 R 基的化学性质可分为<b>极性（亲水）</b>与<b>非极性（疏水）</b>，还有酸性、碱性、含硫、芳香族等。',
      'R 基的极性分布会影响肽链后来的<b>盘曲折叠方式</b>（疏水的 R 基倾向于包在分子内部）。',
      '有 8 种氨基酸人体细胞不能合成，必须从食物中获取，称为<b>必需氨基酸</b>（婴儿为 9 种）。',
    ],
    exam: '高考常考：① 写出氨基酸结构通式；② 根据 R 基判断氨基酸种类；③ 必需氨基酸的概念；④ 人工合成蛋白质时氨基酸的"排列顺序"信息由什么决定（核酸）。',
    build(ctx) {
      const g = new T.Group(); ctx.root.add(g);
      let cur = 'Ala';
      let mol = null, L = [], stat = '';
      function render() {
        if (mol) g.remove(mol.root);
        L.forEach((l) => l.remove()); L = [];
        const data = MOL.buildPeptideAtoms({ n: 1, sidechains: [cur], termini: true });
        centerData(data);
        mol = MOL.renderMol(data, { scale: 1.4 });
        g.add(mol.root);
        const R = data.residues[0];
        const c0 = new V3(0, 0, 0);
        const def = MOL.SIDECHAINS[cur];
        L.push(ctx.lab('主链部分<br><span class="s">（所有氨基酸相同）</span>', outward(data.atoms[R.CA].pos, c0, 3.6), 'c-ink'));
        L.push(ctx.lab('R 基<br><span class="s">' + def.r + '</span>', outward(data.atoms[R.CB].pos, c0, 1.9), 'c-teal'));
        const col = (AA_TYPES.find((a) => a[0] === cur) || [])[1] || '#5eb0a8';
        stat = '<div class="aabox"><div class="aaname">' + def.cn +
          ' <span class="abbr">' + def.abbr + '</span></div>' +
          '<div class="aagrid"><span>R 基</span><b>' + def.r + '</b>' +
          '<span>性质</span><b style="color:' + col + '">' + def.type + '</b></div></div>';
        ctx.stat(stat);
      }
      render();
      ctx.setView(new V3(0, 0, 0), 13, new V3(0.2, 0.25, 1));
      return {
        controls: [{
          type: 'select', label: '选择氨基酸', value: cur,
          options: MOL.DEFAULT_R.map((k) => ({ v: k, t: MOL.SIDECHAINS[k].cn + '（' + k + '）' })),
          onChange(v) { cur = v; render(); },
        }],
        dispose() { L.forEach((l) => l.remove()); },
      };
    },
  };

  /* ============================================================
   * ③ 脱水缩合
   * ============================================================ */
  const s3 = {
    id: 'condense',
    nav: '③ 脱水缩合',
    title: '脱水缩合 —— 形成肽键',
    lead: '一个氨基酸分子的氨基（—NH₂）与另一个氨基酸分子的羧基（—COOH）脱水缩合，脱去一分子水，形成肽键（—CO—NH—）。',
    points: [
      '脱去的水分子中：<b>氧和氢来自羧基（—OH）</b>，另一个<b>氢来自氨基（—H）</b>。',
      '连接两个氨基酸分子的化学键叫<b>肽键</b>，写作 <b>—CO—NH—</b>。',
      '两个氨基酸缩合形成<b>二肽</b>；三个形成三肽；多个形成<b>多肽</b>。',
      '注意：<b>参与肽键形成的只能是 α-碳原子上连接的氨基和羧基</b>，R 基中的氨基、羧基一般不参与。',
    ],
    exam: '必考计算与判断：① n 个氨基酸形成 1 条肽链 → 肽键数 = 脱去水分子数 = n−1；② 形成 m 条肽链 → 肽键数 = n−m；③ 蛋白质相对分子质量 = n·a − 18(n−m)（a 为氨基酸平均相对分子质量）。',
    build(ctx) {
      const g = new T.Group(); ctx.root.add(g);
      // 目标二肽构象
      const full = MOL.buildPeptideAtoms({ n: 2, sidechains: ['Gly', 'Ala'], phi: -120, psi: 120, termini: true });
      // 两个游离氨基酸
      const aaA = MOL.buildPeptideAtoms({ n: 1, sidechains: ['Gly'], termini: true });
      const aaB = MOL.buildPeptideAtoms({ n: 1, sidechains: ['Ala'], termini: true });
      // 把 B 刚体叠合到二肽中第二个残基的位置
      const fit = MOL.rigidFit(
        [aaB.Ns[0].clone(), aaB.CAs[0].clone(), aaB.Cs[0].clone()],
        [full.Ns[1].clone(), full.CAs[1].clone(), full.Cs[1].clone()]
      );
      aaB.atoms.forEach((a) => { a.pos.applyQuaternion(fit.q).add(fit.t); });

      // 重排 B 的氨基两个 H：一个成为成键后的酰胺 H（理想位置），另一个准备脱离
      const NB = aaB.Ns[0], OA = full.Os[0], CA_ = full.Cs[0];
      aaB.atoms[aaB.residues[0].HN[0]].pos = MOL.placeAtom(OA, CA_, NB, 1.01, 119, 180);
      aaB.atoms[aaB.residues[0].HN[1]].pos = MOL.placeAtom(OA, CA_, NB, 1.01, 119, 62);

      // 整体居中
      const all = aaA.atoms.concat(aaB.atoms).map((a) => a.pos);
      const ctr = MOL.centerPoints(all);
      [aaA, aaB].forEach((d) => {
        d.atoms.forEach((a) => a.pos.sub(ctr));
      });
      const N1 = aaB.Ns[0].clone(), C0 = aaA.Cs[0].clone();
      const axis = N1.clone().sub(C0).normalize();

      // 标记脱离的原子
      const RA = aaA.residues[0], RB = aaB.residues[0];
      const iOXT = RA.OXT, iHXT = RA.HXT, iHb = RB.HN[1];
      aaA.atoms[iOXT].detached = true;
      aaA.atoms[iHXT].detached = true;
      aaB.atoms[iHb].detached = true;

      const molA = MOL.renderMol(aaA, { scale: 1.3 });
      const molB = MOL.renderMol(aaB, { scale: 1.3 });
      g.add(molA.root); g.add(molB.root);
      g.attach(molB.atomMeshes[iHb]);

      // 隐藏静态键，改用动态键跟随
      const hideBond = (data, ia, ib) => {
        data.bonds.forEach((b) => {
          if ((b.a === ia && b.b === ib) || (b.a === ib && b.b === ia)) if (b.mesh) b.mesh.visible = false;
        });
      };
      hideBond(aaA, RA.C, iOXT);
      hideBond(aaA, iOXT, iHXT);
      hideBond(aaB, RB.N, iHb);

      const goldMat = new T.MeshPhongMaterial({ color: 0xf0a63c, shininess: 70, specular: 0x664400 });
      const waterMatO = new T.MeshPhongMaterial({ color: MOL.CPK.O, shininess: 60, specular: 0x333333, transparent: true, opacity: 1 });
      const waterMatH = new T.MeshPhongMaterial({ color: MOL.CPK.H, shininess: 60, specular: 0x333333, transparent: true, opacity: 1 });
      const waterBondMat = new T.MeshPhongMaterial({ color: 0xdfe5ec, transparent: true, opacity: 1 });

      const mOXT = molA.atomMeshes[iOXT], mHXT = molA.atomMeshes[iHXT], mHb = molB.atomMeshes[iHb];
      mOXT.material = waterMatO; mHXT.material = waterMatH; mHb.material = waterMatH;

      const dbCO = MOL.dynamicBond(g, molA.atomMeshes[RA.C], mOXT, 0.15);
      const dbOH = MOL.dynamicBond(g, mOXT, mHXT, 0.15);
      const dbNH = MOL.dynamicBond(g, molB.atomMeshes[RB.N], mHb, 0.15);
      const dbW = MOL.dynamicBond(g, mOXT, mHb, 0.15);
      dbW.mesh.material = waterBondMat;
      const dbPep = MOL.dynamicBond(g, molA.atomMeshes[RA.C], molB.atomMeshes[RB.N], 0.17, 0xd9a441);
      dbPep.mesh.material = goldMat;
      dbPep.mesh.visible = false;

      // 高亮光晕
      const glowMatC = new T.MeshBasicMaterial({ color: 0xe8703f, transparent: true, opacity: 0.2 });
      const glowMatN = new T.MeshBasicMaterial({ color: 0x3f7fd0, transparent: true, opacity: 0.2 });
      const glowC = new T.Mesh(new T.SphereGeometry(1.5, 20, 14), glowMatC);
      const glowN = new T.Mesh(new T.SphereGeometry(1.5, 20, 14), glowMatN);
      glowC.position.copy(molA.atomMeshes[RA.C].position);
      glowN.position.copy(N1);
      g.add(glowC); g.add(glowN);

      // 水的目标位置
      const mid = C0.clone().add(N1).multiplyScalar(0.5);
      const upV = new V3(0, 1, 0).addScaledVector(axis, -axis.y).normalize();
      const wPos = mid.clone().addScaledVector(upV, 3.0);
      const hA = new V3(0.78, 0.60, 0), hBv = new V3(-0.78, 0.60, 0);

      const home = {
        OXT: mOXT.position.clone(), HXT: mHXT.position.clone(), Hb: mHb.position.clone(),
      };

      // 标签（原子已居中，方向以原点为基准）
      const O0 = new V3(0, 0, 0);
      const L = [];
      L.push(ctx.lab('羧基 —COOH', outward(aaA.atoms[RA.C].pos, O0, 2.6), 'c-red'));
      const labN = ctx.lab('氨基 —NH₂', outward(N1.clone(), O0, 2.6), 'c-blue');
      L.push(labN);
      const labOH = ctx.lab('—OH 脱下', wPos.clone(), 'c-red');
      const labH = ctx.lab('—H 脱下', wPos.clone(), 'c-blue');
      const labW = ctx.lab('<b>脱去 1 分子 H₂O</b>', wPos.clone(), 'c-water');
      const labPep = ctx.lab('<b>肽键 —CO—NH—</b>', mid.clone().addScaledVector(upV, -1.6), 'c-gold');
      const labDip = ctx.lab('二肽', mid.clone().addScaledVector(upV, 5.2), 'c-ink');
      L.push(labOH, labH, labW, labPep, labDip);
      labOH.hide(); labH.hide(); labW.hide(); labPep.hide(); labDip.hide();

      ctx.setView(new V3(0, 0, 0), 17, new V3(0.1, 0.12, 1));

      let t = 0, playing = true, T0 = 0;
      const statTxt = [
        '<div class="kv"><span>反应前</span><b>2 个游离氨基酸</b></div>',
        '<div class="kv"><span>靠近</span><b>羧基与氨基相互接近</b></div>',
        '<div class="kv"><span>脱水</span><b>脱去 1 分子 H₂O</b></div>',
        '<div class="kv"><span>成键</span><b>形成 1 个肽键 → 二肽</b></div>',
      ];

      function apply(tt) {
        t = Math.min(1, Math.max(0, tt));
        const appr = smoothstep(0.05, 0.72, t);
        const disp = 4.9 * (1 - appr);
        molB.root.position.copy(axis).multiplyScalar(disp);

        const lv = smoothstep(0.24, 0.62, t);
        const float = smoothstep(0.60, 1.0, t) * 3.4;
        mOXT.position.copy(home.OXT).lerp(wPos, lv).addScaledVector(upV, float);
        mHXT.position.copy(home.HXT).lerp(wPos.clone().add(hA), lv).addScaledVector(upV, float);
        mHb.position.copy(home.Hb).addScaledVector(axis, disp)
          .lerp(wPos.clone().add(hBv), lv).addScaledVector(upV, float);

        const fade = 1 - smoothstep(0.78, 1.0, t);
        waterMatO.opacity = fade; waterMatH.opacity = fade; waterBondMat.opacity = fade;

        // 断键
        const broken = t > 0.30;
        dbCO.visible = !broken;
        dbNH.visible = !broken;
        dbOH.visible = true;
        dbW.visible = t > 0.50;

        // 肽键
        const pg = smoothstep(0.70, 0.82, t);
        dbPep.visible = pg > 0.01;
        dbPep.radius = 0.17 * pg;

        // 光晕
        const gp = 1 - smoothstep(0.06, 0.26, t);
        glowMatC.opacity = 0.22 * gp * (0.65 + 0.35 * Math.sin(T0 * 5));
        glowMatN.opacity = 0.22 * gp * (0.65 + 0.35 * Math.sin(T0 * 5));
        glowC.visible = glowN.visible = gp > 0.01;
        glowN.position.copy(N1).addScaledVector(axis, disp);

        g.updateMatrixWorld(true);
        dbCO.update(); dbOH.update(); dbNH.update(); dbW.update(); dbPep.update();

        // 标签
        labN.pos.copy(N1).addScaledVector(axis, disp).addScaledVector(upV, 2.4);
        labOH.show(t > 0.26 && t < 0.60); labOH.pos.copy(mOXT.position).addScaledVector(upV, 1.3);
        labH.show(t > 0.26 && t < 0.60); labH.pos.copy(mHb.position).add(new V3(-0.6, 1.2, 0));
        labW.show(t > 0.48); labW.pos.copy(mOXT.position).addScaledVector(upV, 1.9);
        labPep.show(t > 0.72); labPep.pos.copy(mid).addScaledVector(upV, -2.2);
        labDip.show(t > 0.85);

        ctx.stat('<div class="kvlist">' +
          statTxt[0] + statTxt[1] + statTxt[2] + statTxt[3] +
          '</div><div class="note">拖动下方进度条可以逐步讲解：靠近 → 脱水 → 成键</div>');
      }
      apply(0);
      return {
        controls: [
          { type: 'btn', label: '播放 / 暂停', primary: true, onClick(c) { playing = !playing; } },
          { type: 'btn', label: '重播', onClick(c) { t = 0; playing = true; apply(0); ctrl.range.value = 0; } },
          { type: 'range', min: 0, max: 1, step: 0.005, value: 0, onInput(v) { playing = false; apply(v); } },
        ],
        onControls(cs) { ctrl.range = cs[2].input; ctrl.btn = cs[0]; },
        update(dt, time) {
          T0 = time;
          if (playing) {
            t += dt / 7.2;
            if (t >= 1) { t = 1; playing = false; }
            apply(t);
            if (ctrl.range) ctrl.range.value = t;
          }
        },
        dispose() { L.forEach((l) => l.remove()); },
      };
    },
  };
  // ctrl 句柄（播放控制）
  const ctrl = {};

  /* ============================================================
   * ④ 多肽链
   * ============================================================ */
  const SEQ = ['Ala', 'Ser', 'Gly', 'Leu', 'Glu', 'Val', 'Lys', 'Phe', 'Cys', 'Ala',
    'Gly', 'Ser', 'Val', 'Leu', 'Glu', 'Lys', 'Phe', 'Cys', 'Ala', 'Gly'];
  const PAD_MAIN = ['Gly', 'Ala', 'Cys'];
  const PAD_MORE = ['Ser', 'Val', 'Leu', 'Lys', 'Glu', 'Phe'];
  const s4 = {
    id: 'chain',
    nav: '④ 试一试：搭肽链',
    title: '试一试 —— 用氨基酸搭一条肽链',
    lead: '像课本"试一试"那样：自选氨基酸的种类和数目，让它们一个个脱水缩合连接起来，观察形成的多肽链。',
    points: [
      '每连接一个氨基酸，就<b>脱去 1 分子水、形成 1 个肽键</b>。',
      '<b>肽键数 = 氨基酸数 − 肽链数</b>；脱去的水分子数 = 肽键数。',
      '一端是游离的<b>氨基（N 端）</b>，另一端是游离的<b>羧基（C 端）</b>。',
      '氨基酸的<b>种类、数目、排列顺序</b>可以千变万化 → 蛋白质结构极其多样。',
    ],
    exam: '① 计算：肽键数、脱去水分子数、至少含有的氨基/羧基数、蛋白质相对分子质量；② 至少含有的游离氨基/羧基数 = 肽链数（R 基中的另算）；③ 多肽与蛋白质的区别。',
    build(ctx) {
      const g = new T.Group(); ctx.root.add(g);
      let seq = ['Gly', 'Ala', 'Cys'];
      let mol = null, L = [];
      const extraNH2 = () => seq.filter((k) => k === 'Lys').length;
      const extraCOOH = () => seq.filter((k) => k === 'Glu').length;

      function render() {
        if (mol) { g.remove(mol.root); mol = null; }
        L.forEach((l) => l.remove()); L = [];
        const n = seq.length;
        const data = MOL.buildPeptideAtoms({
          n, phi: -120, psi: 120, sidechains: seq, termini: true, hAlpha: false,
        });
        data.residues.forEach((r) => {
          const col = (AA_TYPES.find((a) => a[0] === r.scKey) || ['', '#5eb0a8'])[1];
          if (r.CB !== undefined) data.atoms[r.CB].color = parseInt(col.slice(1), 16);
        });
        centerData(data);
        mol = MOL.renderMol(data, { scale: 1.2 });
        g.add(mol.root);
        const R0 = data.residues[0], Rn = data.residues[n - 1];
        const c0 = new V3(0, 0, 0);
        L.push(ctx.lab('氨基端（N 端）', outward(data.atoms[R0.N].pos, c0, 2.3), 'c-blue'));
        L.push(ctx.lab('羧基端（C 端）', outward(data.atoms[Rn.C].pos, c0, 2.7), 'c-red'));
        const mid = Math.floor(n / 2);
        L.push(ctx.lab('肽键 —CO—NH—', outward(data.atoms[data.residues[mid].C].pos, c0, 2.5), 'c-gold'));
        L.push(ctx.lab('R 基各不相同', outward(data.atoms[data.residues[Math.min(n - 1, mid + 1)].CB].pos, c0, 2.0), 'c-teal'));
        data.bonds.forEach((b) => {
          if (b.kind === 'peptide' && b.a === data.residues[mid].C && b.mesh) {
            b.mesh.children.forEach((ch) => { ch.material = MOL.mat(0xf0a63c); });
          }
        });
        ctx.stat('<div class="kvlist">' +
          '<div class="kv"><span>氨基酸数</span><b>' + n + '</b></div>' +
          '<div class="kv"><span>肽键数</span><b>' + (n - 1) + '</b></div>' +
          '<div class="kv"><span>脱去水分子数</span><b>' + (n - 1) + '</b></div>' +
          '<div class="kv"><span>至少含游离氨基</span><b>' + (1 + extraNH2()) + '</b></div>' +
          '<div class="kv"><span>至少含游离羧基</span><b>' + (1 + extraCOOH()) + '</b></div>' +
          '</div><div class="note">肽键数 = 氨基酸数 − 肽链数 = 脱去水分子数</div>');
        view();
      }
      function view() {
        ctx.setView(new V3(0, 0, 0), 9 + seq.length * 1.75, new V3(0.05, 0.18, 1));
      }
      render();

      return {
        controls: [{
          type: 'custom',
          mount(box) {
            box.innerHTML = '';
            const mk = (keys, title) => {
              const lb = document.createElement('label');
              lb.className = 'ctllb';
              lb.innerHTML = title;
              box.appendChild(lb);
              const grid = document.createElement('div');
              grid.className = 'aagrid2';
              box.appendChild(grid);
              keys.forEach((k) => {
                const def = MOL.SIDECHAINS[k];
                const col = (AA_TYPES.find((a) => a[0] === k) || ['', '#5eb0a8'])[1];
                const b = document.createElement('button');
                b.className = 'aab';
                b.style.borderColor = col;
                b.innerHTML = '<i style="background:' + col + '"></i><b>' + def.cn +
                  '</b><em>' + k + '</em><s>R = ' + def.r + '</s>';
                b.onclick = () => { if (seq.length < 24) { seq.push(k); render(); drawSeq(); } };
                grid.appendChild(b);
              });
            };
            mk(PAD_MAIN, '点击添加氨基酸（课本"试一试"的三种）');
            mk(PAD_MORE, '更多氨基酸');
            const seqBox = document.createElement('div');
            seqBox.className = 'seqbox';
            box.appendChild(seqBox);
            const acts = document.createElement('div');
            acts.className = 'actrow';
            box.appendChild(acts);
            const mkBtn = (t, fn, cls) => {
              const b = document.createElement('button');
              b.className = 'btn' + (cls ? ' ' + cls : '');
              b.innerHTML = t;
              b.onclick = fn;
              acts.appendChild(b);
              return b;
            };
            mkBtn('删除最后一个', () => { if (seq.length > 1) { seq.pop(); render(); drawSeq(); } });
            mkBtn('随机三个', () => {
              const pool = PAD_MAIN.concat(PAD_MORE);
              seq = [0, 1, 2].map(() => pool[Math.floor(Math.random() * pool.length)]);
              render(); drawSeq();
            });
            mkBtn('清空', () => { seq = ['Gly']; render(); drawSeq(); });
            function drawSeq() {
              seqBox.innerHTML = '';
              seq.forEach((k, i) => {
                const col = (AA_TYPES.find((a) => a[0] === k) || ['', '#5eb0a8'])[1];
                const d = document.createElement('span');
                d.className = 'chip';
                d.style.background = col;
                d.innerHTML = MOL.SIDECHAINS[k].abbr;
                seqBox.appendChild(d);
                if (i < seq.length - 1) {
                  const t = document.createElement('i');
                  t.className = 'tk';
                  t.innerHTML = '';
                  seqBox.appendChild(t);
                }
              });
            }
            drawSeq();
          },
        }],
        dispose() { L.forEach((l) => l.remove()); },
      };
    },
  };

  /* ============================================================
   * ⑤ 二级结构
   * ============================================================ */
  const s5 = {
    id: 'secondary',
    nav: '⑤ 二级结构',
    title: '二级结构 —— α-螺旋与 β-折叠',
    lead: '肽链主链上的羰基氧（C=O）与亚氨基氢（N—H）之间形成氢键，使肽链盘绕成 α-螺旋或折叠成 β-折叠。',
    points: [
      '<b>α-螺旋</b>：肽链盘绕成右手螺旋，<b>每圈约 3.6 个氨基酸残基，螺距 0.54 nm</b>，氢键在<b>第 i 位与第 i+4 位</b>之间。',
      '<b>β-折叠</b>：肽链伸展成锯齿状，若干肽段<b>平行或反平行</b>排列，靠链间氢键连接成片层。',
      '维持二级结构的作用力是<b>氢键</b>，只与<b>主链</b>有关，与 R 基无关。',
      'R 基伸向螺旋外侧或折叠片层的上下两侧，不参与二级结构的形成。',
    ],
    exam: '① 二级结构的类型与氢键位置（α-螺旋为 i 与 i+4）；② 维持各级结构的作用力：二级结构靠氢键，三级结构靠氢键、疏水作用、二硫键等；③ 二级结构只涉及主链构象，不涉及 R 基。',
    build(ctx) {
      const g = new T.Group(); ctx.root.add(g);
      let mode = 'a', mol = null, L = [], hbGroup = null;
      function clear() {
        if (mol) g.remove(mol.root);
        if (hbGroup) g.remove(hbGroup);
        L.forEach((l) => l.remove()); L = [];
      }
      function render() {
        clear();
        if (mode === 'a') {
          const n = 18;
          const data = MOL.buildPeptideAtoms({ n, phi: -57, psi: -47, sidechains: ['Ala'], hAlpha: false });
          centerData(data);
          // 把螺旋轴转到竖直方向，便于观察侧面
          const axisV = data.atoms[data.residues[n - 1].CA].pos.clone()
            .sub(data.atoms[data.residues[0].CA].pos).normalize();
          const q = new T.Quaternion().setFromUnitVectors(axisV, new V3(0, 1, 0));
          data.atoms.forEach((a) => a.pos.applyQuaternion(q));
          mol = MOL.renderMol(data, { scale: 1.1 });
          const hb = findHBonds(data, 3);
          hbGroup = new T.Group();
          drawHBonds(hbGroup, data, hb, 0x1f8a70);
          g.add(hbGroup);
          const c0 = new V3(0, 0, 0);
          L.push(ctx.lab('<b>α-螺旋</b>（右手螺旋）', outward(data.atoms[data.residues[2].CA].pos, c0, 4.2), 'c-ink'));
          L.push(ctx.lab('氢键：第 i 位 C=O ··· 第 i+4 位 N—H',
            outward(data.atoms[data.residues[8].CA].pos, c0, 3.2), 'c-teal'));
          L.push(ctx.lab('R 基伸向螺旋外侧', outward(data.atoms[data.residues[13].CB].pos, c0, 2.8), 'c-purple'));
          L.push(ctx.lab('每圈约 3.6 个残基<br><span class="s">螺距 0.54 nm</span>',
            outward(data.atoms[data.residues[16].CA].pos, c0, 3.6), 'c-gray'));
          ctx.setView(new V3(0, 0, 0), 31, new V3(0.12, 0.08, 1));
          ctx.stat('<div class="kvlist">' +
            '<div class="kv"><span>结构</span><b>α-螺旋</b></div>' +
            '<div class="kv"><span>氢键</span><b>' + hb.length + ' 个（i 与 i+4）</b></div>' +
            '<div class="kv"><span>数据来源</span><b>理想螺旋（键长键角严格）</b></div></div>');
        } else {
          const seg = global.PROTEIN_PDB.segments.hairpin_1gb1;
          const data = buildSegment(seg.residues, {});
          const ctr = MOL.centerPoints(data.atoms.map((a) => a.pos));
          data.atoms.forEach((a) => a.pos.sub(ctr));
          mol = MOL.renderMol(data, { scale: 1.1 });
          const hb = findHBonds(data, 5);
          hbGroup = new T.Group();
          drawHBonds(hbGroup, data, hb, 0x1f8a70);
          g.add(hbGroup);
          const c0 = new V3(0, 0, 0);
          L.push(ctx.lab('<b>β-股 1</b>（第 42–46 位）', outward(data.atoms[data.residues[1].CA].pos, c0, 3.0), 'c-ink'));
          L.push(ctx.lab('<b>β-股 2</b>（第 51–55 位，反向）', outward(data.atoms[data.residues[12].CA].pos, c0, 3.0), 'c-ink'));
          L.push(ctx.lab('链间氢键', outward(data.atoms[data.residues[11].H].pos, c0, 2.6), 'c-teal'));
          L.push(ctx.lab('β-转角', outward(data.atoms[data.residues[7].CA].pos, c0, 2.8), 'c-purple'));
          ctx.setView(new V3(0, 0, 0), 26, new V3(0.1, 0.35, 1));
          ctx.stat('<div class="kvlist">' +
            '<div class="kv"><span>结构</span><b>反平行 β-发夹</b></div>' +
            '<div class="kv"><span>氢键</span><b>' + hb.length + ' 个（股间）</b></div>' +
            '<div class="kv"><span>数据来源</span><b>PDB 1GB1 真实结构</b></div></div>');
        }
        g.add(mol.root);
      }
      render();
      return {
        controls: [{
          type: 'seg', value: 'a',
          options: [{ v: 'a', t: 'α-螺旋' }, { v: 'b', t: 'β-折叠' }],
          onChange(v) { mode = v; render(); },
        }],
        dispose() { L.forEach((l) => l.remove()); },
      };
    },
  };

  /* ============================================================
   * ⑥ 三级结构
   * ============================================================ */
  const s6 = {
    id: 'tertiary',
    nav: '⑥ 三级结构',
    title: '三级结构 —— 肽链盘曲折叠成空间结构',
    lead: '整条肽链在二级结构基础上进一步盘曲、折叠，形成特定的空间结构。只有形成特定空间结构的蛋白质才具有功能。',
    points: [
      '肽链通过<b>盘曲、折叠</b>形成特定的空间结构，主要依靠<b>氢键、疏水作用、二硫键</b>等维持。',
      '<b>蛋白质的空间结构取决于氨基酸序列</b>；序列中一个关键位点的氨基酸改变，功能就可能改变。',
      '有些蛋白质只由<b>一条肽链</b>构成，形成三级结构后即具有功能，如细胞色素 c、肌红蛋白。',
      '实例：<b>肌红蛋白</b> —— 1 条肽链、153 个氨基酸、8 段 α-螺旋，中间结合一个<b>血红素辅基</b>。',
    ],
    exam: '① 结构多样性原因：氨基酸的种类、数目、排列顺序不同，以及肽链盘曲折叠方式及其形成的空间结构不同；② 镰刀型细胞贫血症：血红蛋白 β 链上一个氨基酸被替换 → 空间结构改变 → 功能异常；③ 结构与功能相适应。',
    build(ctx) {
      const g = new T.Group(); ctx.root.add(g);
      const P = global.PROTEIN_PDB.myoglobin;
      const ch = P.chains[0];
      const pts = chainPoints(ch);
      const ctr = MOL.centerPoints(pts);
      pts.forEach((p) => p.sub(ctr));
      const cart = new MOL.Cartoon(pts, { ss: ch.ss, radius: 0.92 });
      g.add(cart.group);

      // 伸展肽链（折叠前的样子）：CA 间距 3.4 Å 的锯齿状长链
      const N = pts.length;
      const ext = [];
      for (let i = 0; i < N; i++) {
        ext.push(new V3((i - (N - 1) / 2) * 3.4, (i % 2 ? 1 : -1) * 1.1, Math.sin(i * 0.55) * 1.6));
      }
      const ec = MOL.centerPoints(ext);
      ext.forEach((p) => p.sub(ec));
      const tmp = pts.map((p) => p.clone());

      const het = P.hets[0];
      const hg = new T.Group();
      MOL.renderHet(hg, het.atoms, { scale: 1.15 });
      const hc = MOL.centerPoints(het.atoms.map((a) => new V3(a.xyz[0], a.xyz[1], a.xyz[2]))).sub(ctr);
      hg.position.copy(ctr).negate();
      g.add(hg);

      const L = [];
      const ss = ch.ss;
      const nHel = (ss.match(/H+/g) || []).length;
      const hi = ss.indexOf('HHHH');
      const ci = ss.indexOf('CCCC');
      const labHel = ctx.lab('<b>α-螺旋</b>（共 ' + nHel + ' 段）', pts[hi + 2].clone(), 'c-orange');
      const labCoil = ctx.lab('无规卷曲 / 回折', pts[ci + 2].clone(), 'c-gray');
      const labHeme = ctx.lab('<b>血红素辅基</b><br><span class="s">含 Fe²⁺，可结合 O₂</span>', hc.clone(), 'c-red');
      const labChain = ctx.lab('1 条肽链 · 153 个氨基酸', pts[2].clone(), 'c-ink');
      const labFold = ctx.lab('<b>肽链盘曲折叠</b><br><span class="s">形成特定空间结构</span>', new V3(0, 14, 0), 'c-teal');
      L.push(labHel, labCoil, labHeme, labChain, labFold);

      ctx.setView(new V3(0, 0, 0), 72, new V3(0.4, 0.25, 1));

      let f = 1, target = 1, playing = false, showHeme = true;
      const handle = { range: null };
      function apply(v) {
        f = Math.min(1, Math.max(0, v));
        const e = f * f * (3 - 2 * f);
        for (let i = 0; i < N; i++) tmp[i].copy(ext[i]).lerp(pts[i], e);
        cart.setPoints(tmp);
        const folded = f > 0.72;
        labHel.show(folded); labCoil.show(folded); labChain.show(folded);
        labFold.show(!folded);
        labHeme.show(showHeme && f > 0.9);
        hg.visible = showHeme && f > 0.9;
        if (handle.range) handle.range.value = f;
        ctx.stat('<div class="kvlist">' +
          '<div class="kv"><span>实例</span><b>肌红蛋白（PDB 1MBO）</b></div>' +
          '<div class="kv"><span>肽链数</span><b>1 条 · 153 个氨基酸</b></div>' +
          '<div class="kv"><span>当前状态</span><b>' +
          (f < 0.05 ? '伸展的肽链（一级结构）' : f > 0.95 ? '特定空间结构（有功能）' : '正在盘曲折叠…') +
          '</b></div>' +
          '<div class="kv"><span>二级结构</span><b>' + (folded ? nHel + ' 段 α-螺旋 + 回折' : '尚未形成') + '</b></div>' +
          '</div><div class="note">橙红色 = α-螺旋，灰色 = 无规卷曲。拖动进度条可逐步讲解折叠过程</div>');
      }
      apply(1);
      return {
        controls: [
          { type: 'btn', label: '播放折叠 / 复原', primary: true, onClick() { target = f > 0.5 ? 0 : 1; playing = true; } },
          {
            type: 'range', min: 0, max: 1, step: 0.005, value: 1,
            onInput(v) { playing = false; apply(v); },
          },
          { type: 'toggle', label: '显示血红素辅基', value: true, onChange(v) { showHeme = v; apply(f); } },
        ],
        onControls(cs) { handle.range = cs[1].input; },
        update(dt) {
          if (!playing) return;
          const sp = dt / 4.2;
          f += Math.sign(target - f) * Math.min(sp, Math.abs(target - f));
          apply(f);
          if (Math.abs(target - f) < 0.001) playing = false;
        },
        dispose() { L.forEach((l) => l.remove()); cart.dispose(); },
      };
    },
  };

  /* ============================================================
   * ⑦ 四级结构
   * ============================================================ */
  const s7 = {
    id: 'quaternary',
    nav: '⑦ 四级结构',
    title: '四级结构 —— 多条肽链聚集成功能蛋白',
    lead: '有些蛋白质由多条肽链组成：每条肽链各有自己的三级结构，再通过一定的方式聚集在一起，形成具有功能的蛋白质。',
    points: [
      '实例：<b>血红蛋白</b>由 <b>4 条肽链（2 条 α 链 + 2 条 β 链）</b>和 <b>4 个血红素</b>组成。',
      '每条肽链结合一个血红素，每个血红素中央有一个 <b>Fe²⁺</b>，可结合 1 分子 O₂。',
      '因此<b>1 个血红蛋白分子最多可运输 4 分子 O₂</b>。',
      '只有四条肽链都齐全、结构完整时，血红蛋白才能正常运输氧气 —— 体现<b>结构与功能相适应</b>。',
    ],
    exam: '① 血红蛋白的组成（4 条肽链 + 4 个血红素）与运输 O₂ 的数量关系；② 四级结构的概念（多条肽链聚集）；③ 结合糖尿病、镰刀型贫血等实例分析"结构改变 → 功能改变"。',
    build(ctx) {
      const g = new T.Group(); ctx.root.add(g);
      const P = global.PROTEIN_PDB.hemoglobin;
      const CHAIN_COL = [0xe8703f, 0x3f7fd0, 0xe8703f, 0x3f7fd0];
      const all = [];
      P.chains.forEach((ch) => { ch.xyz.forEach((p) => all.push(new V3(p[0], p[1], p[2]))); });
      const ctr = MOL.centerPoints(all);
      const carts = [], centers = [];
      P.chains.forEach((ch, i) => {
        const pts = chainPoints(ch);
        pts.forEach((p) => p.sub(ctr));
        const cart = new MOL.Cartoon(pts, { color: CHAIN_COL[i], radius: 0.92 });
        carts.push(cart);
        g.add(cart.group);
        centers.push(MOL.centerPoints(pts));
      });
      const hetG = new T.Group(); g.add(hetG);
      P.hets.forEach((h) => {
        const sub = new T.Group();
        MOL.renderHet(sub, h.atoms, { scale: 1.1 });
        MOL.centerPoints(h.atoms.map((a) => new V3(a.xyz[0], a.xyz[1], a.xyz[2]))).sub(ctr);
        sub.position.copy(ctr).negate();
        hetG.add(sub);
      });
      const o2g = new T.Group(); g.add(o2g); o2g.visible = false;
      P.hets.forEach((h) => {
        const fe = h.atoms.find((a) => a.el === 'Fe');
        if (!fe) return;
        const p = new V3(fe.xyz[0], fe.xyz[1], fe.xyz[2]).sub(ctr);
        const d = p.clone().normalize().multiplyScalar(2.4);
        const a = p.clone().add(d).addScaledVector(new V3(0, 1, 0), 0.7);
        const b = p.clone().add(d).addScaledVector(new V3(0, 1, 0), -0.7).addScaledVector(d.clone().normalize(), 0.6);
        MOL.addSphere(o2g, 'O', a, 0.42, { color: 0xe04b3a });
        MOL.addSphere(o2g, 'O', b, 0.42, { color: 0xe04b3a });
        MOL.cylBetween(o2g, a, b, 0.16, 0xe04b3a);
      });

      const L = [];
      L.push(ctx.lab('<b>α 链</b>（141 个氨基酸）×2', centers[0].clone(), 'c-orange'));
      L.push(ctx.lab('<b>β 链</b>（146 个氨基酸）×2', centers[1].clone(), 'c-blue'));
      L.push(ctx.lab('<b>血红素</b>（含 Fe²⁺）×4', centers[0].clone().multiplyScalar(0.35).add(new V3(0, 9, 0)), 'c-red'));
      ctx.setView(new V3(0, 0, 0), 56, new V3(0.35, 0.3, 1));
      ctx.stat('<div class="kvlist">' +
        '<div class="kv"><span>实例</span><b>血红蛋白（PDB 2HHB）</b></div>' +
        '<div class="kv"><span>肽链</span><b>4 条（2α + 2β）</b></div>' +
        '<div class="kv"><span>氨基酸总数</span><b>' + P.chains.reduce((s, c) => s + c.n, 0) + '</b></div>' +
        '<div class="kv"><span>辅基</span><b>4 个血红素</b></div>' +
        '<div class="kv"><span>最多运输</span><b>4 分子 O₂</b></div></div>');
      return {
        controls: [{
          type: 'toggle', label: '显示结合的 O₂（示意）', value: false,
          onChange(v) { o2g.visible = v; },
        }],
        dispose() { L.forEach((l) => l.remove()); carts.forEach((c) => c.dispose()); },
      };
    },
  };

  /* ============================================================
   * ⑧ 二硫键（胰岛素）
   * ============================================================ */
  const s9 = {
    id: 'disulfide',
    nav: '⑧ 二硫键',
    title: '二硫键 —— 把空间结构"锁"住',
    lead: '两个半胱氨酸（Cys）R 基上的巯基（—SH）各脱去一个 H，形成二硫键（—S—S—），把肽链的不同部位连接固定，稳定蛋白质的空间结构。',
    points: [
      '二硫键形成于两个<b>半胱氨酸</b>的 R 基（—CH₂SH）之间，写作 <b>—S—S—</b>，形成时脱去 2 个 H。',
      '<b>胰岛素</b>由 <b>2 条肽链</b>组成：A 链 21 个氨基酸、B 链 30 个氨基酸，两条链之间靠 <b>2 个链间二硫键</b>连接（A7–B7、A20–B19）。',
      'A 链内部还有 <b>1 个链内二硫键</b>（A6–A11），共 3 个二硫键。',
      '二硫键是维持<b>空间结构</b>的重要作用力；二硫键被破坏 → 空间结构改变 → <b>功能丧失</b>。',
    ],
    exam: '① 二硫键的形成部位与写法（—S—S—，两个 Cys 的 R 基之间）；② 胰岛素两条链：肽键数 = 氨基酸数 − 肽链数 = 51 − 2 = 49；③ 破坏二硫键只破坏空间结构，肽键未断裂，仍能与双缩脲试剂产生紫色反应。',
    build(ctx) {
      const g = new T.Group(); ctx.root.add(g);
      const P = global.PROTEIN_PDB.insulin;
      const CH_COL = [0xe8703f, 0x3f7fd0];
      const all = [];
      P.chains.forEach((ch) => { ch.xyz.forEach((p) => all.push(new V3(p[0], p[1], p[2]))); });
      const ctr = MOL.centerPoints(all);

      const carts = [], centers = [], ptsList = [];
      P.chains.forEach((ch, i) => {
        const pts = chainPoints(ch);
        pts.forEach((p) => p.sub(ctr));
        ptsList.push(pts);
        const cart = new MOL.Cartoon(pts, { color: CH_COL[i], radius: 0.82 });
        carts.push(cart);
        g.add(cart.group);
        centers.push(MOL.centerPoints(pts));
      });
      // B 链被"推开"的方向
      const sepDir = centers[1].clone().sub(centers[0]);
      sepDir.y = Math.abs(sepDir.y) + 4;
      sepDir.normalize();

      // 二硫键：硫原子随各自肽链移动，中间的键单独一组
      const bondG = new T.Group(); g.add(bondG);
      const SSCOL = 0xf2c14e, SSCOL2 = 0xe0a92c;
      const info = [];
      P.ssbonds.forEach((b) => {
        const pa = new V3(b.ca[0], b.ca[1], b.ca[2]).sub(ctr);
        const pb = new V3(b.cb[0], b.cb[1], b.cb[2]).sub(ctr);
        const inter = b.a[0] !== b.b[0];
        const gi = P.chains.findIndex((c) => c.id === b.a[0]);
        const gj = P.chains.findIndex((c) => c.id === b.b[0]);
        const gp = carts[gi] ? carts[gi].group : g;
        const gq = carts[gj] ? carts[gj].group : g;
        const sa = MOL.addSphere(gp, 'S', pa, 0.46, { color: SSCOL });
        const sb = MOL.addSphere(gq, 'S', pb, 0.46, { color: SSCOL });
        // 把硫原子接回对应残基的 CA，避免悬浮在主链之外。
        const sideLinks = [];
        [[b.a, pa, gp, b.betaA], [b.b, pb, gq, b.betaB]].forEach(([id, sulfur, parent, beta]) => {
          const chain = P.chains.find((c) => c.id === id[0]);
          const idx = Number(id.slice(1)) - chain.first;
          if (idx >= 0 && beta) {
            const cb = new V3(...beta).sub(ctr);
            sideLinks.push(MOL.cylBetween(parent, new V3(...chain.xyz[idx]).sub(ctr), cb, 0.12, SSCOL));
            sideLinks.push(MOL.cylBetween(parent, cb, sulfur, 0.12, SSCOL));
          }
        });
        const link = MOL.dynamicBond(bondG, sa, sb, inter ? 0.17 : 0.14, inter ? SSCOL : SSCOL2);
        info.push({ pa: pa.clone(), pb: pb.clone(), gi, gj, sa, sb, sideLinks, inter,
          name: b.a + '–' + b.b, d: b.d, link });
      });

      const L = [];
      L.push(ctx.lab('<b>A 链</b>（21 个氨基酸）', centers[0].clone().add(new V3(0, 4.4, 0)), 'c-orange'));
      const labB = ctx.lab('<b>B 链</b>（30 个氨基酸）', centers[1].clone().add(new V3(0, -4.6, 0)), 'c-blue');
      L.push(labB);
      const ssLabels = [];
      info.forEach((b) => {
        const mid = b.pa.clone().add(b.pb).multiplyScalar(0.5);
        const t = '<b>' + (b.inter ? '链间' : '链内') + '二硫键</b> —S—S—<br><span class="s">' +
          b.name.replace('A', 'A 链 ').replace('B', 'B 链 ') + ' · ' + b.d + ' Å</span>';
        ssLabels.push(ctx.lab(t, mid, b.inter ? 'c-gold' : 'c-purple'));
      });
      L.push.apply(L, ssLabels);
      const labBreak = ctx.lab('<b>二硫键断裂</b><br><span class="s">两条链分开 → 空间结构破坏 → 功能丧失</span>',
        new V3(0, 12, 0), 'c-red');
      labBreak.hide();
      L.push(labBreak);

      ctx.setView(new V3(0, 0, 0), 52, new V3(0.3, 0.25, 1));

      let sep = 0, target = 0, showSS = true;
      function stat() {
        ctx.stat('<div class="kvlist">' +
          '<div class="kv"><span>实例</span><b>胰岛素（PDB 4INS）</b></div>' +
          '<div class="kv"><span>肽链</span><b>2 条（A 链 21 + B 链 30）</b></div>' +
          '<div class="kv"><span>氨基酸总数</span><b>51</b></div>' +
          '<div class="kv"><span>肽键数</span><b>51 − 2 = 49</b></div>' +
          '<div class="kv"><span>二硫键</span><b>' + (sep > 0 ? '0（已断裂）' : '3 个（链间 2 + 链内 1）') + '</b></div>' +
          '<div class="kv"><span>状态</span><b>' + (sep > 0 ? '<span style="color:#d0402c">空间结构破坏 · 功能丧失</span>' : '空间结构稳定 · 功能正常') + '</b></div>' +
          '</div><div class="note">黄色球 = 半胱氨酸 R 基上的硫原子（S），黄色棒 = 二硫键 —S—S—</div>');
      }
      function apply(v) {
        sep = Math.min(1, Math.max(0, v));
        const e = sep * sep * (3 - 2 * sep);
        carts[1].group.position.copy(sepDir).multiplyScalar(e * 16);
        carts[0].group.position.copy(sepDir).multiplyScalar(-e * 3);
        // 所有二硫键同步断开；显示开关只改变可见性。
        g.updateMatrixWorld(true);
        info.forEach((b, i) => {
          b.link.visible = showSS && sep === 0;
          b.link.update();
          b.sa.visible = b.sb.visible = showSS;
          b.sideLinks.forEach(mesh => { mesh.visible = showSS; });
          ssLabels[i].show(showSS && sep === 0);
          ssLabels[i].pos.copy(b.pa).add(carts[b.gi].group.position)
            .add(b.pb).add(carts[b.gj].group.position).multiplyScalar(0.5);
        });
        labBreak.show(sep > 0);
        L[0].pos.copy(centers[0]).add(new V3(0, 4.4, 0)).add(carts[0].group.position);
        labB.pos.copy(centers[1]).add(new V3(0, -4.6, 0)).add(carts[1].group.position);
        stat();
      }
      apply(0);
      return {
        controls: [
          {
            type: 'seg', value: 'intact',
            options: [{ v: 'intact', t: '完整（有二硫键）' }, { v: 'broken', t: '断开二硫键' }],
            onChange(v) { target = (v === 'broken') ? 1 : 0; },
          },
          { type: 'toggle', label: '显示二硫键', value: true, onChange(v) { showSS = v; apply(sep); } },
        ],
        update(dt) {
          if (Math.abs(sep - target) < 0.0004) return;
          const sp = dt / 1.6;
          sep += Math.sign(target - sep) * Math.min(sp, Math.abs(target - sep));
          apply(sep);
        },
        dispose() { L.forEach((l) => l.remove()); carts.forEach((c) => c.dispose()); },
      };
    },
  };

  /* ============================================================
   * ⑨ 结构与功能 / 变性
   * ============================================================ */
  const s8 = {
    id: 'function',
    nav: '⑨ 结构与功能',
    title: '结构决定功能 —— 蛋白质变性与水解',
    lead: '蛋白质的功能与其特定的空间结构紧密相关。空间结构被破坏，功能就会丧失。',
    points: [
      '<b>结构多样性</b>：氨基酸的种类、数目、排列顺序不同，加上肽链盘曲折叠方式及形成的空间结构不同 → 蛋白质结构极其多样。',
      '<b>功能多样性</b>：结构蛋白（肌肉、毛发）、催化（酶）、运输（血红蛋白）、调节（胰岛素）、免疫（抗体）。',
      '<b>变性</b>：高温、强酸、强碱、重金属盐等破坏蛋白质的<b>空间结构</b>，但<b>肽键没有断裂</b>（一级结构不变），功能丧失，一般不可逆。',
      '<b>水解</b>：在酸、碱或酶的作用下<b>肽键断裂</b>，蛋白质被分解为多肽、氨基酸，一级结构被破坏。',
    ],
    exam: '高频考点：① 变性的本质（空间结构破坏、肽键不断、不可逆）与实例（煮鸡蛋、高温灭菌、酒精消毒、重金属中毒喝牛奶）；② 变性后仍能与双缩脲试剂产生紫色反应（肽键还在！）；③ 盐析（可逆）与变性（不可逆）的区别。',
    build(ctx) {
      const g = new T.Group(); ctx.root.add(g);
      const P = global.PROTEIN_PDB.myoglobin;
      const ch = P.chains[0];
      const pts = chainPoints(ch);
      const ctr = MOL.centerPoints(pts);
      pts.forEach((p) => p.sub(ctr));
      const coil = MOL.randomCoilCA(pts.length, 20261013);
      const c2 = MOL.centerPoints(coil);
      coil.forEach((p) => p.sub(c2));
      // 与折叠态的质心对齐
      const c3 = MOL.centerPoints(coil);
      coil.forEach((p) => p.sub(c3));

      const cart = new MOL.Cartoon(pts, { ss: ch.ss, radius: 0.92 });
      g.add(cart.group);
      const frags = new T.Group(); g.add(frags); frags.visible = false;
      const fragCarts = [];
      for (let k = 0; k < 4; k++) {
        const a = Math.floor((pts.length * k) / 4), b = Math.floor((pts.length * (k + 1)) / 4);
        const seg = pts.slice(a, b);
        const fc = new MOL.Cartoon(seg, { ss: ch.ss.slice(a, b), radius: 0.92 });
        const off = new V3((k - 1.5) * 11, (k % 2 ? 3 : -3), 0);
        fc.group.position.copy(off);
        frags.add(fc.group);
        fragCarts.push(fc);
      }

      const het = P.hets[0];
      const hg = new T.Group();
      MOL.renderHet(hg, het.atoms, { scale: 1.1 });
      const hc = MOL.centerPoints(het.atoms.map((a) => new V3(a.xyz[0], a.xyz[1], a.xyz[2]))).sub(ctr);
      hg.position.copy(ctr).negate();
      g.add(hg);

      const L = [];
      const lbl = ctx.lab('<b>天然构象</b><br><span class="s">具有特定空间结构 → 有功能</span>', new V3(0, -13, 0), 'c-teal');
      const lbl2 = ctx.lab('<b>变性</b>：空间结构被破坏，肽键未断裂', new V3(0, -15, 0), 'c-red');
      const lbl3 = ctx.lab('<b>水解</b>：肽键断裂 → 一级结构被破坏', new V3(0, -15, 0), 'c-red');
      L.push(lbl, lbl2, lbl3);
      lbl2.hide(); lbl3.hide();

      ctx.setView(new V3(0, 0, 0), 74, new V3(0.3, 0.2, 1));

      let mode = 'native', m = 0, target = 0;
      const tmp = [];
      for (let i = 0; i < pts.length; i++) tmp.push(new V3());
      function stat() {
        const map = {
          native: ['<div class="kv"><span>状态</span><b>天然构象</b></div><div class="kv"><span>空间结构</span><b>完整</b></div><div class="kv"><span>肽键</span><b>完好</b></div><div class="kv"><span>功能</span><b>正常</b></div>', ''],
          denature: ['<div class="kv"><span>状态</span><b>变性</b></div><div class="kv"><span>空间结构</span><b style="color:#d0402c">被破坏</b></div><div class="kv"><span>肽键</span><b>完好（一级结构不变）</b></div><div class="kv"><span>功能</span><b style="color:#d0402c">丧失</b></div>', '<div class="note">高温、强酸、强碱、重金属盐等 → 空间结构改变，一般不可逆</div>'],
          hydrolyze: ['<div class="kv"><span>状态</span><b>水解</b></div><div class="kv"><span>空间结构</span><b style="color:#d0402c">被破坏</b></div><div class="kv"><span>肽键</span><b style="color:#d0402c">断裂</b></div><div class="kv"><span>产物</span><b>多肽 → 氨基酸</b></div>', '<div class="note">酸、碱或蛋白酶催化下肽键断裂，一级结构被破坏</div>'],
        };
        ctx.stat('<div class="kvlist">' + map[mode][0] + '</div>' + map[mode][1]);
      }
      stat();
      return {
        controls: [{
          type: 'seg', value: 'native',
          options: [{ v: 'native', t: '天然构象' }, { v: 'denature', t: '变性' }, { v: 'hydrolyze', t: '水解' }],
          onChange(v) {
            mode = v;
            target = (v === 'denature') ? 1 : 0;
            frags.visible = (v === 'hydrolyze');
            cart.group.visible = (v !== 'hydrolyze');
            hg.visible = (v !== 'hydrolyze');
            lbl.show(v === 'native'); lbl2.show(v === 'denature'); lbl3.show(v === 'hydrolyze');
            stat();
          },
        }],
        update(dt) {
          if (Math.abs(m - target) > 0.0005) {
            const spd = dt / 1.6;
            m += Math.sign(target - m) * Math.min(spd, Math.abs(target - m));
            const e = m * m * (3 - 2 * m);
            for (let i = 0; i < pts.length; i++) tmp[i].copy(pts[i]).lerp(coil[i], e);
            cart.setPoints(tmp);
          }
        },
        dispose() { L.forEach((l) => l.remove()); cart.dispose(); fragCarts.forEach((f) => f.dispose()); },
      };
    },
  };

  global.SCENES = [s1, s2, s3, s4, s5, s6, s7, s9, s8];
})(window);
