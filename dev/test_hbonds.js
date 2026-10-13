// 检验氢键自动识别：理想 α-螺旋 / 真实 α-螺旋(1GB1) / 真实反平行 β-发夹(1GB1)
global.window = global;
global.THREE = require('../vendor/three.min.js');
require('../js/mol.js');
const MOL = global.MOL;
const V3 = THREE.Vector3;
require('../data/proteins.js');
const PDB = global.PROTEIN_PDB;

function hbonds(HN, O, minSep) {
  const out = [];
  for (const h of HN) for (const o of O) {
    if (Math.abs(h.res - o.res) < minSep) continue;
    const d = h.pos.distanceTo(o.pos);
    if (d < 2.5) out.push({ h: h.res, o: o.res, d: +d.toFixed(2) });
  }
  return out;
}

// ---------- 1. 理想 α-螺旋 ----------
{
  const st = MOL.buildPeptideAtoms({ n: 14, phi: -57, psi: -47, sidechains: false, showH: true });
  const HN = [], O = [];
  st.residues.forEach((r, i) => {
    (r.HN || []).forEach((h) => HN.push({ res: i, pos: st.atoms[h].pos }));
    O.push({ res: i, pos: st.atoms[r.O].pos });
  });
  const hb = hbonds(HN, O, 3);
  console.log(`理想 α-螺旋 (14残基): 氢键 ${hb.length} 个`);
  console.log('  ', hb.map((b) => `H${b.h}···O${b.o}(${b.d})`).join(' '));
}

// ---------- 2. 真实结构片段 ----------
function testSegment(key, label, minSep) {
  const seg = PDB.segments[key];
  const res = seg.residues;
  const O = res.map((r, i) => ({ res: i, pos: new V3(r.O[0], r.O[1], r.O[2]) }));
  const HN = [];
  res.forEach((r, i) => {
    if (i === 0) return;
    const p = res[i - 1];
    const h = MOL.placeAtom(
      new V3(p.O[0], p.O[1], p.O[2]),
      new V3(p.C[0], p.C[1], p.C[2]),
      new V3(r.N[0], r.N[1], r.N[2]),
      1.01, 119, 180
    );
    HN.push({ res: i, pos: h });
  });
  const hb = hbonds(HN, O, minSep);
  console.log(`${label} (${res.length}残基, ${seg.note}): 氢键 ${hb.length} 个`);
  console.log('  ', hb.map((b) => `H${b.h}···O${b.o}(${b.d})`).join(' '));
  // 相邻残基 CA 距离，确认主链连续
  let mx = 0;
  for (let i = 1; i < res.length; i++) {
    const a = new V3(res[i - 1].CA[0], res[i - 1].CA[1], res[i - 1].CA[2]);
    const b = new V3(res[i].CA[0], res[i].CA[1], res[i].CA[2]);
    mx = Math.max(mx, a.distanceTo(b));
  }
  console.log('   相邻 CA 最大间距:', mx.toFixed(2), '(应≈3.8，说明主链连续)');
}
testSegment('helix_1gb1', '真实 α-螺旋', 3);
console.log();
testSegment('hairpin_1gb1', '真实反平行 β-发夹', 3);
