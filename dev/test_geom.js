// 验证 NeRF（自然延伸参考系）生成肽主链与 α-螺旋的正确性
const THREE = require('../vendor/three.min.js');
const V = THREE.Vector3;
const rad = (d) => (d * Math.PI) / 180;

function placeAtom(A, B, C, r, thetaDeg, phiDeg) {
  const bcDir = C.clone().sub(B).normalize();
  const abDir = B.clone().sub(A);
  let n = new V().crossVectors(bcDir, abDir);
  if (n.lengthSq() < 1e-8) n.set(0, 0, 1).cross(bcDir);
  n.normalize();
  const m = new V().crossVectors(n, bcDir).normalize();
  const t = rad(thetaDeg), p = rad(phiDeg);
  const e1 = bcDir.clone().negate();
  return new V()
    .addScaledVector(e1, r * Math.cos(t))
    .addScaledVector(m, -r * Math.cos(p) * Math.sin(t))
    .addScaledVector(n, r * Math.sin(p) * Math.sin(t))
    .add(C);
}

function dist(a, b) { return a.distanceTo(b); }
function angle(a, b, c) {
  const v1 = a.clone().sub(b), v2 = c.clone().sub(b);
  return (Math.acos(v1.dot(v2) / (v1.length() * v2.length())) * 180) / Math.PI;
}
function dihedral(a, b, c, d) {
  const b1 = b.clone().sub(a), b2 = c.clone().sub(b), b3 = d.clone().sub(c);
  const n1 = new V().crossVectors(b1, b2).normalize();
  const n2 = new V().crossVectors(b2, b3).normalize();
  const m1 = new V().crossVectors(n1, b2.clone().normalize());
  const x = n1.dot(n2), y = m1.dot(n2);
  return (Math.atan2(y, x) * 180) / Math.PI;
}

// --- 基本正确性 ---
const A = new V(1, 0, 0), B = new V(0, 0, 0), C = new V(0.5, 0.9, 0);
const D = placeAtom(A, B, C, 1.5, 110, 75);
console.log('键长 1.5 ->', dist(C, D).toFixed(4));
console.log('键角 110 ->', angle(B, C, D).toFixed(3));
console.log('二面角 75 ->', dihedral(A, B, C, D).toFixed(3));

// --- 构建肽主链 ---
function buildBackbone(n, phi, psi) {
  const atoms = [];
  // 种子 N1, CA1, C1
  const N1 = new V(0, 0, 0);
  const CA1 = new V(1.458, 0, 0);
  const C1 = new V(1.458 + 1.525 * Math.cos(rad(180 - 111)), 1.525 * Math.sin(rad(180 - 111)), 0);
  atoms.push({ N: N1, CA: CA1, C: C1 });
  let prev = atoms[0];
  for (let i = 1; i < n; i++) {
    const N = placeAtom(prev.N, prev.CA, prev.C, 1.329, 116, psi);
    const CA = placeAtom(prev.CA, prev.C, N, 1.458, 121, 180);
    const C = placeAtom(prev.C, N, CA, 1.525, 111, phi);
    const cur = { N, CA, C };
    atoms.push(cur);
    prev = cur;
  }
  return atoms;
}

// --- α-螺旋: phi=-57, psi=-47 ---
const hel = buildBackbone(20, -57, -47);
console.log('\n=== α-螺旋 ===');
console.log('肽键 C-N 长度:', dist(hel[0].C, hel[1].N).toFixed(3), '(应≈1.33)');
console.log('α-碳间距 CA(i)-CA(i+1):', dist(hel[0].CA, hel[1].CA).toFixed(3), '(应≈3.8)');
// 每圈残基数 & 螺距: 用 i 与 i+4 的 CA 距离判断
console.log('CA(i)-CA(i+4):', dist(hel[3].CA, hel[7].CA).toFixed(3), '(α-螺旋应≈6.2)');
console.log('CA(i)-CA(i+3):', dist(hel[3].CA, hel[6].CA).toFixed(3), '(应≈5.0-5.5)');
// 手性判据: 右手 α-螺旋的 CA 虚拟二面角约为 +50°
function caTorsion(atoms, i) {
  return dihedral(atoms[i - 1].CA, atoms[i].CA, atoms[i + 1].CA, atoms[i + 2].CA);
}
let tor = 0, cnt = 0;
for (let i = 5; i < 15; i++) { tor += caTorsion(hel, i); cnt++; }
console.log('CA 虚拟二面角(右手α螺旋≈+50°):', (tor / cnt).toFixed(2));
const helL = buildBackbone(20, 57, 47);
let torL = 0;
for (let i = 5; i < 15; i++) torL += dihedral(helL[i - 1].CA, helL[i].CA, helL[i + 1].CA, helL[i + 2].CA);
console.log('  镜像对照 (phi=+57,psi=+47):', (torL / cnt).toFixed(2));

// 每圈残基数 / 螺距: 用 5~15 号残基拟合螺旋轴
function helixParams(atoms, lo, hi) {
  const pts = atoms.slice(lo, hi + 1).map((a) => a.CA);
  const cen = new V();
  pts.forEach((p) => cen.add(p));
  cen.divideScalar(pts.length);
  // 主轴 = 首尾连线
  const ax = pts[pts.length - 1].clone().sub(pts[0]).normalize();
  let ang = 0, n = 0;
  const ref = pts[1].clone().sub(cen); ref.addScaledVector(ax, -ref.dot(ax)).normalize();
  let prevAng = 0, unwrapped = 0, prev = ref.clone();
  for (let i = 2; i < pts.length; i++) {
    const p = pts[i].clone().sub(cen); p.addScaledVector(ax, -p.dot(ax)).normalize();
    const cross = new V().crossVectors(prev, p);
    let d = Math.acos(Math.max(-1, Math.min(1, prev.dot(p)))) * 180 / Math.PI;
    if (cross.dot(ax) < 0) d = -d;
    unwrapped += d; prev = p; n++;
  }
  ang = Math.abs(unwrapped) / n;
  const rise = pts[pts.length - 1].clone().sub(pts[0]).dot(ax) / (pts.length - 1);
  return { perRes: ang, resPerTurn: 360 / ang, rise, pitch: (360 / ang) * Math.abs(rise) };
}
const hp = helixParams(hel, 4, 16);
console.log('每残基转角:', hp.perRes.toFixed(1), '° (应≈100°)');
console.log('每圈残基数:', hp.resPerTurn.toFixed(2), '(应≈3.6)');
console.log('每残基上升:', hp.rise.toFixed(3), 'Å (应≈1.5)');
console.log('螺距:', hp.pitch.toFixed(2), 'Å (应≈5.4)');

// --- 伸展链: phi=-120, psi=+120 ---
const ext = buildBackbone(12, -120, 120);
console.log('\n=== β-伸展链 ===');
console.log('CA(i)-CA(i+1):', dist(ext[0].CA, ext[1].CA).toFixed(3), '(应≈3.8)');
console.log('CA(i)-CA(i+2):', dist(ext[3].CA, ext[5].CA).toFixed(3), '(应≈6.7)');
console.log('首尾直线距离(12残基):', dist(ext[0].CA, ext[11].CA).toFixed(2), '(应≈36)');
let maxdev = 0;
const line = ext[11].CA.clone().sub(ext[0].CA).normalize();
for (const a of ext) {
  const v = a.CA.clone().sub(ext[0].CA); v.addScaledVector(line, -v.dot(line));
  maxdev = Math.max(maxdev, v.length());
}
console.log('偏离主轴最大值:', maxdev.toFixed(2), '(锯齿状应≈1)');
