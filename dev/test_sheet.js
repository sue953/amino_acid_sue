// 反平行 β-折叠：绕片层法线旋转 180° 生成第二条股，扫描股间距与轴向错位
global.window = global;
global.THREE = require('../vendor/three.min.js');
require('../js/mol.js');
const MOL = global.MOL;
const V3 = THREE.Vector3;

const st = MOL.buildPeptideAtoms({ n: 10, phi: -120, psi: 120, sidechains: false, showH: true });
const CA = st.CAs;
const HN = [], O = [];
st.residues.forEach((r, i) => {
  (r.HN || []).forEach((h) => HN.push({ res: i, pos: st.atoms[h].pos }));
  O.push({ res: i, pos: st.atoms[r.O].pos });
});

// 局部坐标系
const u = CA[CA.length - 1].clone().sub(CA[0]).normalize();
const w = new V3().crossVectors(CA[1].clone().sub(CA[0]), CA[2].clone().sub(CA[1])).normalize();
const v = new V3().crossVectors(w, u).normalize();
console.log('股轴 u =', u.toArray().map((x) => x.toFixed(2)).join(','));
console.log('法线 w =', w.toArray().map((x) => x.toFixed(2)).join(','));
console.log('褶皱 v =', v.toArray().map((x) => x.toFixed(2)).join(','));

const R = new THREE.Matrix4().makeRotationAxis(w, Math.PI);
const T = (p, du, dv) => p.clone().applyMatrix4(R).addScaledVector(u, du).addScaledVector(v, dv);

const results = [];
for (let du = -3.0; du <= 3.0; du += 0.2) {
  for (let dv = 4.0; dv <= 6.0; dv += 0.1) {
    const HNb = HN.map((h) => ({ res: h.res, pos: T(h.pos, du, dv) }));
    const Ob = O.map((o) => ({ res: o.res, pos: T(o.pos, du, dv) }));
    let n = 0, sum = 0;
    const pairs = [];
    for (const h of HN) for (const o of Ob) {
      const d = h.pos.distanceTo(o.pos);
      if (d < 2.35) { n++; sum += d; pairs.push(`A:H${h.res}-B:O${o.res}=${d.toFixed(2)}`); }
    }
    for (const h of HNb) for (const o of O) {
      const d = h.pos.distanceTo(o.pos);
      if (d < 2.35) { n++; sum += d; pairs.push(`B:H${h.res}-A:O${o.res}=${d.toFixed(2)}`); }
    }
    if (n >= 5) results.push({ du: +du.toFixed(1), dv: +dv.toFixed(1), n, avg: +(sum / n).toFixed(2), pairs });
  }
}
results.sort((a, b) => b.n - a.n);
console.log('\n候选（>=5 对）:');
results.slice(0, 8).forEach((r) => {
  console.log(`  du=${r.du} dv=${r.dv}  对数=${r.n} 平均=${r.avg}`);
  console.log('    ', r.pairs.slice(0, 8).join('  '));
});
if (!results.length) console.log('  无');
