(function () {
  'use strict';
  function contrast(gap, h, imbalance, threshold) {
    const p2 = 0.10, q2 = 0.10, p1 = p2 + gap, q1 = q2 + imbalance;
    const low1 = p1 - h * q1, low2 = p2 - h * q2;
    const valid = [low1, low1 + h, low2, low2 + h, q1].every(x => x >= -1e-12 && x <= 1 + 1e-12);
    return { valid, corrected: gap - h * imbalance, threshold,
      tipping: h > 0 ? gap / h : null,
      decisionTipping: h > 0 ? (gap - threshold) / h : null };
  }
  if (typeof module !== 'undefined' && module.exports) module.exports = { contrast };
  if (typeof document === 'undefined') return;
  function setup() {
    const root = document.getElementById('a10-tipping-simulator');
    if (!root) return;
    const field = name => root.querySelector('[name="' + name + '"]');
    const svg = root.querySelector('svg'), result = root.querySelector('[data-result]');
    const ns = 'http://www.w3.org/2000/svg';
    function el(tag, attrs, content) {
      const n = document.createElementNS(ns, tag);
      Object.entries(attrs).forEach(([k, v]) => n.setAttribute(k, v));
      if (content !== undefined) n.textContent = content;
      svg.appendChild(n); return n;
    }
    function update() {
      const gap = +field('gap').value, h = +field('h').value,
        imbalance = +field('imbalance').value, threshold = +field('threshold').value;
      root.querySelectorAll('input').forEach(n => {
        root.querySelector('[data-value="' + n.name + '"]').textContent = (100 * +n.value).toFixed(1) + ' points';
      });
      const r = contrast(gap, h, imbalance, threshold);
      root._a10Result = r;
      result.textContent = r.valid ? 'Standardized difference: ' + (100 * r.corrected).toFixed(1) +
        ' percentage points (' + (1000 * r.corrected).toFixed(0) + ' additional expected events per 1,000). ' +
        (r.corrected >= threshold - 1e-12 ? 'Meets' : 'Falls below') + ' the chosen decision threshold. ' +
        (r.tipping === null ? 'With h = 0, changing latent prevalence does not change the contrast.' :
          'The algebraic sign tipping point is an imbalance of ' + (100 * r.tipping).toFixed(1) + ' points; assess whether it is plausible.') :
        'This combination implies an underlying event probability outside [0,1]. Reduce the latent effect or imbalance; no valid model result is reported.';
      svg.replaceChildren();
      const x = v => 65 + v / .6 * 570, y = v => 260 - (v + .25) / .37 * 215;
      el('title', {}, 'Standardized probability difference across latent-prevalence imbalances');
      el('line', {x1:65,y1:y(0),x2:635,y2:y(0),stroke:'#666','stroke-dasharray':'4 4'});
      el('line', {x1:65,y1:y(threshold),x2:635,y2:y(threshold),stroke:'#b45309','stroke-dasharray':'7 3'});
      for (const v of [0,.2,.4,.6]) {
        el('text',{x:x(v),y:283,'text-anchor':'middle','font-size':13},String(100*v));
      }
      for (const v of [-.2,-.1,0,.1]) {
        el('text',{x:53,y:y(v)+4,'text-anchor':'end','font-size':13},String(100*v));
      }
      let path = '';
      for (let i=0;i<=120;i++) {
        const imbalanceAt = i / 200, c = contrast(gap,h,imbalanceAt,threshold);
        if (c.valid) path += (path ? ' L ' : 'M ') + x(imbalanceAt) + ' ' + y(c.corrected);
      }
      el('path',{d:path,fill:'none',stroke:'#2563eb','stroke-width':3});
      if (r.valid) el('circle',{cx:x(imbalance),cy:y(r.corrected),r:6,fill:'#15803d'});
      el('text',{x:350,y:315,'text-anchor':'middle','font-size':14},'Latent-prevalence imbalance (percentage points)');
      el('text',{x:65,y:22,'font-size':14},'Standardized difference (percentage points); orange = decision threshold');
      root.dataset.ready = 'true';
    }
    root.addEventListener('input',update);
    root.querySelector('button').addEventListener('click',() => {
      const defaults={gap:.06,h:.2,imbalance:.2,threshold:.02};
      Object.entries(defaults).forEach(([k,v])=>field(k).value=v);update();
    });
    update();
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded',setup);
  else setup();
})();
