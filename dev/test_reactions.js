// Exercise the actual scene callbacks and Three.js geometry without a GPU.
const assert = require('node:assert/strict');
global.window = global;
global.THREE = require('../vendor/three.min.js');
require('../js/mol.js');
require('../data/proteins.js');
require('../js/scenes.js');
function build(id) {
  const labels = [];
  const ctx = {
    root: new THREE.Group(), html: '',
    lab(text, pos) {
      const item = { text, pos: pos.clone(), visible: true,
        show(v) { this.visible = v !== false; }, hide() { this.visible = false; }, remove() {} };
      labels.push(item); return item;
    },
    stat(html) { this.html = html; }, setView() {},
  };
  return { ctx, labels, api: SCENES.find(s => s.id === id).build(ctx) };
}
const reactionScene = SCENES.find(s => s.title.includes('脱水缩合'));
const reaction = build(reactionScene.id);
const rg = reaction.ctx.root.children[0];
const a = rg.children[0], b = rg.children[1];
const slider = reaction.api.controls.find(c => c.type === 'range');
slider.onInput(1);
// Backbone spheres retain chemical bond lengths after rigid fit and centering.
for (const mol of [a, b]) {
  const atoms = mol.children[0].children.filter(o => o.isMesh);
  assert.ok(Math.abs(atoms[0].position.distanceTo(atoms[1].position) - 1.458) < 1e-6);
  assert.ok(Math.abs(atoms[1].position.distanceTo(atoms[2].position) - 1.525) < 1e-6);
}
const pos = mesh => mesh.getWorldPosition(new THREE.Vector3());
assert.ok(Math.abs(pos(a.children[0].children[2]).distanceTo(pos(b.children[0].children[0])) - 1.329) < 1e-6);
const links = rg.children.filter(o => o.isMesh && o.geometry.type === 'CylinderGeometry');
assert.equal(links.length, 5);
assert.deepEqual(links.map(o => o.visible), [false, true, false, true, true]);
slider.onInput(0);
assert.deepEqual(links.map(o => o.visible), [true, true, true, false, false]);
const ss = build('disulfide');
const sg = ss.ctx.root.children[0];
const bonds = sg.children[2].children;
assert.equal(bonds.length, 2);
// Four cysteine sidechains each expose CA, CB and two hydrogen spheres.
const sideSpheres = sg.children.slice(0, 2).flatMap(c => c.children)
  .filter(o => o.isMesh && o.geometry.type === 'SphereGeometry');
assert.equal(sideSpheres.length, 20);
assert.ok(sideSpheres.every(o => Number.isFinite(o.position.x)));
assert.equal(ss.labels.filter(l => l.text.includes('链内二硫键</b>')).length, 1);
assert.equal(ss.labels.filter(l => l.text.includes('链间二硫键</b>')).length, 1);
assert.ok(bonds.every(b => b.visible));
ss.api.controls[1].onChange(false);
assert.ok(bonds.every(b => !b.visible));
assert.equal(ss.labels.find(l => l.text.includes('二硫键断裂')).visible, false);
assert.ok(ss.ctx.html.includes('2 个（链内 1 + 链间 1）'));
ss.api.controls[1].onChange(true);
ss.api.controls[0].onChange('broken');
for (let i = 0; i < 40; i++) ss.api.update(0.05);
assert.ok(bonds.every(b => !b.visible));
assert.ok(ss.ctx.html.includes('0（已断裂）'));
ss.api.controls[0].onChange('intact');
for (let i = 0; i < 40; i++) ss.api.update(0.05);
assert.ok(bonds.every(b => b.visible));
assert.equal(ss.labels.find(l => l.text.includes('二硫键断裂')).visible, false);
console.log('Reaction geometry, bond visibility, disulfide toggle and restore checks passed.');
