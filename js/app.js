/* ============================================================
 * app.js —— 页面外壳、场景调度与交互
 * ============================================================ */
(function () {
  'use strict';
  const T = window.THREE;
  const V3 = T.Vector3;
  const SCENES = window.SCENES;

  let renderer, scene, camera, controls, clock;
  let root, viewport, labelsEl;
  let labels = [];
  let api = null, curIdx = -1, tween = null;
  let spin = true, showLabels = true;
  let elapsed = 0;

  /* ---------------- 标签层 ---------------- */
  const ctx = {
    root: null,
    spin: 1,
    lab(text, pos, cls) {
      const el = document.createElement('div');
      el.className = 'lab ' + (cls || '');
      el.innerHTML = text;
      labelsEl.appendChild(el);
      const item = {
        el, pos: pos.clone(),
        show(v) { el.style.display = v === false ? 'none' : ''; },
        hide() { this.show(false); },
        remove() { el.remove(); const i = labels.indexOf(item); if (i >= 0) labels.splice(i, 1); },
      };
      labels.push(item);
      return item;
    },
    stat(html) { document.getElementById('statBox').innerHTML = html; },
    setView(center, dist, dir) {
      const to = center.clone();
      const pos = center.clone().addScaledVector(dir.clone().normalize(), dist);
      tween = {
        fp: camera.position.clone(), tp: pos,
        ft: controls.target.clone(), tt: to, t: 0, dur: 0.75,
      };
    },
  };

  function updateLabels() {
    const w = viewport.clientWidth, h = viewport.clientHeight;
    for (const it of labels) {
      if (!showLabels) { it.el.style.display = 'none'; continue; }
      if (it._hidden) { it.el.style.display = 'none'; continue; }
      it.el.style.display = '';
      const p = it.pos.clone();
      root.localToWorld(p);
      const v = p.project(camera);
      if (v.z > 1) { it.el.style.display = 'none'; continue; }
      const x = (v.x * 0.5 + 0.5) * w, y = (-v.y * 0.5 + 0.5) * h;
      if (x < -200 || x > w + 200 || y < -100 || y > h + 100) { it.el.style.display = 'none'; continue; }
      it.el.style.transform = 'translate(-50%,-50%) translate(' + x.toFixed(1) + 'px,' + y.toFixed(1) + 'px)';
    }
  }
  // show(v) 需要记录状态，避免与全局开关打架
  const _show = (v) => {};
  function patchShow(item) {
    const orig = item.show;
    item.show = function (v) { item._hidden = (v === false); orig(v); };
    return item;
  }

  /* ---------------- 场景切换 ---------------- */
  function goto(i) {
    if (i < 0 || i >= SCENES.length) return;
    const S = SCENES[i];
    // 清理
    if (api && api.dispose) api.dispose();
    while (root.children.length) root.remove(root.children[0]);
    labels.slice().forEach((l) => l.remove());
    root.rotation.set(0, 0, 0);
    document.getElementById('statBox').innerHTML = '';
    document.getElementById('ctlBox').innerHTML = '';

    // 面板
    document.getElementById('iTitle').innerHTML = S.title;
    document.getElementById('iLead').innerHTML = S.lead;
    document.getElementById('iPoints').innerHTML = S.points.map((p) => '<li>' + p + '</li>').join('');
    document.getElementById('iExam').innerHTML = S.exam;

    try {
      api = S.build(ctx) || {};
    } catch (e) {
      api = {};
      window.__buildErr = (S.id || '?') + ': ' + ((e && e.stack) || e);
      console.error(window.__buildErr);
    }
    curIdx = i;

    // 控件
    if (api.controls && api.controls.length) buildControls(api.controls);
    if (api.onControls) api.onControls(api._made || []);

    document.querySelectorAll('.step').forEach((el, k) => el.classList.toggle('on', k === i));
    document.getElementById('btnPrev').disabled = (i === 0);
    document.getElementById('btnNext').disabled = (i === SCENES.length - 1);
    const cur = document.querySelector('.step.on');
    if (cur) cur.scrollIntoView({ block: 'nearest' });
  }

  function buildControls(specs) {
    const box = document.getElementById('ctlBox');
    box.innerHTML = '';
    const made = [];
    specs.forEach((sp) => {
      const row = document.createElement('div');
      row.className = 'ctl ctl-' + sp.type;
      if (sp.label && sp.type !== 'btn') {
        const lb = document.createElement('label');
        lb.className = 'ctllb';
        lb.innerHTML = sp.label;
        row.appendChild(lb);
      }
      const made0 = { spec: sp, el: row };
      if (sp.type === 'btn') {
        const b = document.createElement('button');
        b.className = 'btn' + (sp.primary ? ' primary' : '');
        b.innerHTML = sp.label;
        b.onclick = () => sp.onClick(made0);
        row.appendChild(b);
        made0.button = b;
      } else if (sp.type === 'range') {
        const wrap = document.createElement('div');
        wrap.className = 'rng';
        const inp = document.createElement('input');
        inp.type = 'range'; inp.min = sp.min; inp.max = sp.max; inp.step = sp.step; inp.value = sp.value;
        const val = document.createElement('span');
        val.className = 'rngval'; val.innerHTML = sp.value;
        inp.oninput = () => { val.innerHTML = inp.value; sp.onInput(parseFloat(inp.value)); };
        wrap.appendChild(inp); wrap.appendChild(val);
        row.appendChild(wrap);
        made0.input = inp; made0.val = val;
      } else if (sp.type === 'select') {
        const sel = document.createElement('select');
        sp.options.forEach((o) => {
          const op = document.createElement('option');
          op.value = o.v; op.innerHTML = o.t;
          sel.appendChild(op);
        });
        sel.value = sp.value;
        sel.onchange = () => sp.onChange(sel.value);
        row.appendChild(sel);
        made0.input = sel;
      } else if (sp.type === 'seg') {
        const wrap = document.createElement('div');
        wrap.className = 'seg';
        sp.options.forEach((o) => {
          const b = document.createElement('button');
          b.innerHTML = o.t;
          b.className = o.v === sp.value ? 'on' : '';
          b.onclick = () => {
            wrap.querySelectorAll('button').forEach((x) => x.classList.remove('on'));
            b.classList.add('on');
            sp.onChange(o.v);
          };
          wrap.appendChild(b);
        });
        row.appendChild(wrap);
      } else if (sp.type === 'custom') {
        if (sp.mount) sp.mount(row);
      } else if (sp.type === 'toggle') {
        const lab = document.createElement('label');
        lab.className = 'sw';
        const cb = document.createElement('input');
        cb.type = 'checkbox'; cb.checked = !!sp.value;
        const sp2 = document.createElement('span');
        sp2.className = 'swt';
        lab.appendChild(cb); lab.appendChild(sp2);
        lab.appendChild(document.createTextNode(sp.label));
        cb.onchange = () => sp.onChange(cb.checked);
        row.appendChild(lab);
        made0.input = cb;
      }
      box.appendChild(row);
      made.push(made0);
    });
    api._made = made;
  }

  /* ---------------- 主循环 ---------------- */
  function animate() {
    requestAnimationFrame(animate);
    const dt = Math.min(0.05, clock.getDelta());
    elapsed += dt;
    if (tween) {
      tween.t += dt / tween.dur;
      const k = Math.min(1, tween.t);
      const e = k < 0.5 ? 2 * k * k : 1 - Math.pow(-2 * k + 2, 2) / 2;
      camera.position.lerpVectors(tween.fp, tween.tp, e);
      controls.target.lerpVectors(tween.ft, tween.tt, e);
      camera.lookAt(controls.target);
      if (k >= 1) tween = null;
    } else {
      controls.update();
    }
    if (spin) root.rotation.y += dt * 0.16;
    if (api && api.update) api.update(dt, elapsed);
    updateLabels();
    renderer.render(scene, camera);
  }

  /* ---------------- 初始化 ---------------- */
  function init() {
    viewport = document.getElementById('viewport');
    labelsEl = document.getElementById('labels');

    renderer = new T.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    renderer.setSize(viewport.clientWidth, viewport.clientHeight);
    renderer.outputEncoding = T.sRGBEncoding;
    viewport.insertBefore(renderer.domElement, labelsEl);

    scene = new T.Scene();
    scene.fog = new T.Fog(0xeef1f6, 90, 260);
    camera = new T.PerspectiveCamera(42, viewport.clientWidth / viewport.clientHeight, 0.1, 2000);
    camera.position.set(0, 0, 30);

    controls = new T.OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.08;
    controls.minDistance = 4;
    controls.maxDistance = 320;

    const hemi = new T.HemisphereLight(0xffffff, 0xd3dced, 0.62);
    scene.add(hemi);
    const d1 = new T.DirectionalLight(0xffffff, 0.62);
    d1.position.set(9, 14, 16); scene.add(d1);
    const d2 = new T.DirectionalLight(0xffffff, 0.26);
    d2.position.set(-12, -5, -10); scene.add(d2);
    scene.add(new T.AmbientLight(0xffffff, 0.24));

    root = new T.Group();
    ctx.root = root;
    scene.add(root);
    window.__dbg = { scene, root, camera, T };

    clock = new T.Clock();

    // 步骤导航
    const nav = document.getElementById('steps');
    SCENES.forEach((S, i) => {
      const b = document.createElement('button');
      b.className = 'step';
      b.innerHTML = '<i>' + (i + 1) + '</i><span>' + S.nav.replace(/^[①②③④⑤⑥⑦⑧⑨]\s*/, '') + '</span>';
      b.onclick = () => goto(i);
      nav.appendChild(b);
    });

    document.getElementById('btnPrev').onclick = () => goto(curIdx - 1);
    document.getElementById('btnNext').onclick = () => goto(curIdx + 1);

    const bs = document.getElementById('btnSpin');
    bs.onclick = () => { spin = !spin; bs.classList.toggle('off', !spin); bs.innerHTML = (spin ? '自动旋转：开' : '自动旋转：关'); };
    const bl = document.getElementById('btnLabel');
    bl.onclick = () => {
      showLabels = !showLabels;
      bl.classList.toggle('off', !showLabels);
      bl.innerHTML = (showLabels ? '标签：显示' : '标签：隐藏');
    };
    document.getElementById('btnFull').onclick = () => {
      const el = document.getElementById('stage');
      if (!document.fullscreenElement) el.requestFullscreen && el.requestFullscreen();
      else document.exitFullscreen();
    };
    document.addEventListener('fullscreenchange', onResize);

    window.addEventListener('resize', onResize);
    window.addEventListener('keydown', (e) => {
      if (e.target && /INPUT|SELECT|TEXTAREA/.test(e.target.tagName)) return;
      if (e.key === 'ArrowRight') goto(curIdx + 1);
      else if (e.key === 'ArrowLeft') goto(curIdx - 1);
      else if (e.key === ' ') {
        e.preventDefault();
        if (api && api._made) {
          const b = api._made.find((m) => m.spec.type === 'btn');
          if (b) b.button.click();
        }
      }
    });

    // lab 的 show 增强
    const _lab = ctx.lab;
    ctx.lab = function (t, p, c) { return patchShow(_lab(t, p, c)); };

    goto(0);
    animate();
  }

  function onResize() {
    if (!renderer) return;
    const w = viewport.clientWidth, h = viewport.clientHeight;
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
    renderer.setSize(w, h);
  }

  window.addEventListener('DOMContentLoaded', init);
})();
