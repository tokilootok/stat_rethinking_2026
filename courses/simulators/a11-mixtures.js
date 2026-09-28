/* A11: deterministic distribution explorer. No network, dependencies or fitting. */
(() => {
  'use strict';
  const logistic = x => 1 / (1 + Math.exp(-x));
  const logGamma = z => {
    const c = [676.5203681218851, -1259.1392167224028, 771.3234287776531,
      -176.6150291621406, 12.507343278686905, -.13857109526572012,
      9.984369578019572e-6, 1.5056327351493116e-7];
    if (z < .5) return Math.log(Math.PI) - Math.log(Math.sin(Math.PI*z)) - logGamma(1-z);
    z -= 1; let a = .99999999999980993;
    c.forEach((v, i) => { a += v/(z+i+1); });
    const t = z + 7.5;
    return .5*Math.log(2*Math.PI)+(z+.5)*Math.log(t)-t+Math.log(a);
  };
  const logBeta = (a, b) => logGamma(a)+logGamma(b)-logGamma(a+b);
  const poisson = (y, m) => Math.exp(-m+y*Math.log(m)-logGamma(y+1));
  const models = {
    bb: {
      title: 'Beta-binomial: shared cohort probability',
      controls: [['mu', 'Mean probability μ', .02, .98, .01, .2], ['kappa', 'Concentration κ', .5, 80, .5, 10], ['n', 'Opportunities n', 2, 60, 1, 30]],
      calculate: v => {
        const probabilities = Array.from({length:v.n+1}, (_,y) => Math.exp(logGamma(v.n+1)-logGamma(y+1)-logGamma(v.n-y+1)+logBeta(y+v.mu*v.kappa,v.n-y+(1-v.mu)*v.kappa)-logBeta(v.mu*v.kappa,(1-v.mu)*v.kappa)));
        return {probabilities, mean:v.n*v.mu, variance:v.n*v.mu*(1-v.mu)*(v.n+v.kappa)/(1+v.kappa), note:'Larger κ reduces variation between cohorts; the mean stays fixed when μ and n stay fixed.'};
      }
    },
    nb: {
      title: 'Gamma-Poisson: heterogeneous rates',
      controls: [['mu', 'Mean count μ', .2, 15, .2, 4], ['phi', 'Gamma shape φ', .2, 30, .2, 2]],
      calculate: v => ({probabilities:Array.from({length:61}, (_,y) => Math.exp(logGamma(y+v.phi)-logGamma(v.phi)-logGamma(y+1)+v.phi*Math.log(v.phi/(v.phi+v.mu))+y*Math.log(v.mu/(v.phi+v.mu)))),mean:v.mu,variance:v.mu+v.mu*v.mu/v.phi,note:'Smaller φ adds both a heavier tail and more zeros, without an offline state.'})
    },
    zip: {
      title: 'Zero-inflated Poisson: two routes to zero',
      controls: [['psi', 'Offline probability ψ', 0, .95, .01, .3], ['lambda', 'Active-state mean λ', .1, 10, .1, 2]],
      calculate: v => ({probabilities:Array.from({length:36}, (_,y) => (y===0?v.psi:0)+(1-v.psi)*poisson(y,v.lambda)),mean:(1-v.psi)*v.lambda,variance:(1-v.psi)*v.lambda*(1+v.psi*v.lambda),note:`P(offline | zero) = ${(v.psi/(v.psi+(1-v.psi)*Math.exp(-v.lambda))).toFixed(3)}. A zero does not reveal its source.`})
    },
    ordinal: {
      title: 'Ordered outcome: cumulative probabilities become category probabilities',
      controls: [['eta', 'Linear predictor η', -4, 4, .1, 0], ['center', 'Cutpoint center', -2, 2, .1, 0], ['gap', 'Distance between cutpoints', .2, 5, .1, 2]],
      calculate: v => {
        const f1=logistic(v.center-v.gap/2-v.eta), f2=logistic(v.center+v.gap/2-v.eta);
        return {probabilities:[f1,f2-f1,1-f2],labels:['Low','Medium','High'],note:`Cutpoints ${(v.center-v.gap/2).toFixed(2)} and ${(v.center+v.gap/2).toFixed(2)}. Cumulative probabilities: ${f1.toFixed(3)}, ${f2.toFixed(3)}, 1. Positive η moves probability toward High.`};
      }
    },
    monotonic: {
      title: 'Ordered predictor: allocate an endpoint effect across unequal steps',
      controls: [['beta','Total log-odds shift β',-3,3,.1,1.5],['w1','Relative weight: step 1',.1,5,.1,1],['w2','Relative weight: step 2',.1,6,.1,6],['w3','Relative weight: step 3',.1,5,.1,3]],
      calculate: v => {
        const sum=v.w1+v.w2+v.w3, d=[v.w1/sum,v.w2/sum,v.w3/sum];
        return {probabilities:[0,d[0],d[0]+d[1],1].map(s=>logistic(v.beta*s-1)),labels:['Level 1','Level 2','Level 3','Level 4'],isCurve:true,note:`Normalized increments δ = (${d.map(x=>x.toFixed(3)).join(', ')}). Bars are P(High) at different levels; they do not sum to one. Cutpoints are fixed at −1 and 1.`};
      }
    }
  };
  // Expose pure calculations for numerical verification, including outside a browser.
  if (typeof module !== 'undefined') module.exports = models;
  if (typeof document === 'undefined') return;
  const mount = document.getElementById('a11-explorer');
  if (!mount) return;
  mount.innerHTML = '<label for="a11-family"><strong>Model to explore</strong></label> <select id="a11-family"></select><div id="a11-controls"></div><p id="a11-values" aria-live="polite"></p><div id="a11-plot"></div><p id="a11-note"></p><p><small>These are distributions at chosen parameter values, not fitted posteriors. Executed Bayesian fits appear in the lesson examples.</small></p>';
  const select = mount.querySelector('select');
  Object.entries(models).forEach(([key, model]) => {const o=document.createElement('option');o.value=key;o.textContent=model.title;select.appendChild(o);});
  const render = () => {
    const model = models[select.value], values={};
    model.controls.forEach(([key])=>{values[key]=Number(mount.querySelector('#a11-'+key).value);mount.querySelector('#a11-'+key+'-value').textContent=values[key];});
    const r=model.calculate(values), p=r.probabilities;
    let footer='';
    if (r.mean!==undefined) footer=`Mean = ${r.mean.toFixed(3)}; variance = ${r.variance.toFixed(3)}; P(0) = ${p[0].toFixed(3)}.`;
    else if (!r.isCurve) footer=`Category probabilities: ${p.map(x=>x.toFixed(3)).join(', ')}; sum = ${p.reduce((a,b)=>a+b,0).toFixed(6)}.`;
    else footer=`P(High): ${p.map(x=>x.toFixed(3)).join(', ')}.`;
    if (!r.labels && select.value!=='bb') footer+=` Probability beyond the plotted range = ${Math.max(0,1-p.reduce((a,b)=>a+b,0)).toFixed(4)}.`;
    mount.querySelector('#a11-values').textContent=footer;
    mount.querySelector('#a11-note').textContent=r.note;
    const width=760,height=285,left=45,top=20,base=245,plotWidth=690,step=plotWidth/p.length;
    const ymax=r.isCurve?1:Math.min(1,Math.max(...p)*1.15);
    let svg=`<svg viewBox="0 0 ${width} ${height}" role="img" aria-label="${model.title}" style="width:100%;max-width:900px;background:#fff;color:#111"><title>${model.title}</title>`;
    [0,.5,1].forEach(f=>{const y=base-(base-top)*f;svg+=`<line x1="${left}" y1="${y}" x2="735" y2="${y}" stroke="#d1d5db"/><text x="2" y="${y+4}" font-size="12">${(f*ymax).toFixed(2)}</text>`;});
    p.forEach((prob,i)=>{const h=prob/ymax*(base-top),x=left+i*step;svg+=`<rect x="${x+1}" y="${base-h}" width="${Math.max(1,step-2)}" height="${h}" fill="#246b8e"><title>${r.labels?r.labels[i]:i}: ${prob.toFixed(5)}</title></rect>`;if(r.labels||i%Math.max(1,Math.ceil(p.length/12))===0)svg+=`<text x="${x+step/2}" y="265" text-anchor="middle" font-size="12">${r.labels?r.labels[i]:i}</text>`;});
    mount.querySelector('#a11-plot').innerHTML=svg+'</svg>';
  };
  const setup=()=>{
    const box=mount.querySelector('#a11-controls');box.innerHTML='';
    models[select.value].controls.forEach(([key,label,min,max,step,value])=>{
      const row=document.createElement('div');row.style.cssText='display:flex;align-items:center;gap:1em;flex-wrap:wrap;margin:0.7em 0';
      row.innerHTML=`<label for="a11-${key}" style="min-width:230px">${label}</label><input id="a11-${key}" type="range" min="${min}" max="${max}" step="${step}" value="${value}"><output id="a11-${key}-value" for="a11-${key}">${value}</output>`;
      box.appendChild(row);row.querySelector('input').addEventListener('input',render);
    });render();
  };
  select.addEventListener('change',setup);setup();
})();
