/* Local, dependency-free teaching simulators. Grid quadrature, not MCMC. */
(function () {
  'use strict';
  const sigmoid = x => x >= 0 ? 1 / (1 + Math.exp(-x)) : Math.exp(x) / (1 + Math.exp(x));
  const softplus = x => Math.max(x, 0) + Math.log1p(Math.exp(-Math.abs(x)));
  const normal = (x, m, s) => Math.exp(-0.5 * ((x - m) / s) ** 2) / (s * Math.sqrt(2 * Math.PI));
  const axis = (lo, hi, n) => Array.from({length:n}, (_,i) => lo + (i + 0.5) * (hi - lo) / n);
  const sum = a => a.reduce((s,x) => s+x, 0);
  function normalize(logw) {
    let peak = -Infinity;
    for (const x of logw) peak = Math.max(peak,x);
    const w = logw.map(x => Math.exp(x-peak)), total = sum(w);
    return w.map(x => x/total);
  }
  function stats(x,w) {
    const mean = sum(x.map((v,i) => v*w[i]));
    const variance = sum(x.map((v,i) => (v-mean)**2*w[i]));
    let acc=0, lower=x[0], upper=x[x.length-1], found=false;
    for(let i=0;i<x.length;i++) {
      acc+=w[i];
      if(!found && acc>=0.05) {lower=x[i];found=true;}
      if(acc>=0.95) {upper=x[i];break;}
    }
    return {mean,variance,sd:Math.sqrt(variance),lower,upper};
  }
  function intercept(ma,sa,n,y) {
    const x=axis(Math.min(-12,ma-7*sa),Math.max(12,ma+7*sa),1201);
    const logprior=x.map(a=>-0.5*((a-ma)/sa)**2);
    const prior=normalize(logprior);
    const posterior=normalize(x.map((a,i)=>logprior[i]+y*a-n*softplus(a)));
    const p=x.map(sigmoid);
    return {x,p,prior,posterior,alpha:stats(x,posterior),prob:stats(p,posterior),priorProb:stats(p,prior)};
  }
  function regression(ma,sa,mb,sb,data,size=201,bounds=null) {
    const a=axis(bounds?bounds[0]:Math.min(-8,ma-6*sa),bounds?bounds[1]:Math.max(8,ma+6*sa),size);
    const b=axis(bounds?bounds[2]:Math.min(-8,mb-6*sb),bounds?bounds[3]:Math.max(8,mb+6*sb),size);
    const logs=[];
    for(const av of a) for(const bv of b) {
      let lp=-0.5*((av-ma)/sa)**2-0.5*((bv-mb)/sb)**2;
      for(const d of data) {const eta=av+bv*d.x;lp+=d.y*eta-d.n*softplus(eta);}
      logs.push(lp);
    }
    const w=normalize(logs),wa=a.map(()=>0),wb=b.map(()=>0);
    let edge=0,acc=0,next=0;
    const draws=[],count=2400;
    for(let i=0;i<size;i++) for(let j=0;j<size;j++) {
      const weight=w[i*size+j];wa[i]+=weight;wb[j]+=weight;
      if(i<3 || j<3 || i>=size-3 || j>=size-3) edge+=weight;
      acc+=weight;
      while(next<count && acc>=(next+0.5)/count) {draws.push([a[i],b[j]]);next++;}
    }
    while(draws.length<count) draws.push([a[size-1],b[size-1]]);
    const ast=stats(a,wa),bst=stats(b,wb);
    if(!bounds && ((a[1]-a[0])>ast.sd/3 || (b[1]-b[0])>bst.sd/3)) {
      const aw=8*Math.max(ast.sd,a[1]-a[0]),bw=8*Math.max(bst.sd,b[1]-b[0]);
      return regression(ma,sa,mb,sb,data,size,[ast.mean-aw,ast.mean+aw,bst.mean-bw,bst.mean+bw]);
    }
    let covariance=0;
    for(let i=0;i<size;i++) for(let j=0;j<size;j++) covariance+=w[i*size+j]*(a[i]-ast.mean)*(b[j]-bst.mean);
    const curve=axis(-2.04,2.04,51).map(x=>{
      const pp=draws.map(([av,bv])=>sigmoid(av+bv*x)).sort((u,v)=>u-v);
      const s=Math.hypot(sa,sb*x),m=ma+mb*x;
      return {x,mean:sum(pp)/count,lo:pp[Math.floor(.05*count)],hi:pp[Math.floor(.95*count)],
        priorMedian:sigmoid(m),priorLo:sigmoid(m-1.644854*s),priorHi:sigmoid(m+1.644854*s)};
    });
    return {a,b,wa,wb,alpha:ast,beta:bst,correlation:covariance/(ast.sd*bst.sd),curve,edge};
  }
  function simulate(a,b,n,seed) {
    let state=seed>>>0;
    function random() {state=(Math.imul(1664525,state)+1013904223)>>>0;return state/4294967296;}
    return Array.from({length:9},(_,i)=>{
      const x=-2+i*.5,p=sigmoid(a+b*x);let y=0;
      for(let k=0;k<n;k++) if(random()<p)y++;
      return {x,n,y};
    });
  }
  const api={sigmoid,normal,stats,intercept,regression,simulate};
  if(typeof module!=='undefined' && module.exports) module.exports=api;
  if(typeof document==='undefined') return;

  const blue='#176caa',red='#be3d36',ink='#263445';
  const fmt=(x,d=3)=>Number(x).toFixed(d);
  function field(id,label,value,min,max,step=0.1,type='range') {
    return `<label>${label}<span class="a09-field"><input data-key="${id}" type="${type}" min="${min}" max="${max}" step="${step}" value="${value}" aria-label="${label}"><output>${value}</output></span></label>`;
  }
  function chart(root,title,label) {
    const box=document.createElement('figure');
    box.innerHTML=`<figcaption>${title}</figcaption><canvas width="760" height="290" role="img" aria-label="${label}"></canvas>`;
    root.append(box);return box.querySelector('canvas');
  }
  function draw(canvas,series,{xmin,xmax,ymax=1,xlabel='',ylabel='Densité',points=[]}) {
    const c=canvas.getContext('2d'),W=760,H=290,L=62,R=20,T=12,B=49;
    const X=x=>L+(x-xmin)/(xmax-xmin)*(W-L-R),Y=y=>H-B-y/ymax*(H-T-B);
    c.clearRect(0,0,W,H);c.fillStyle='#fff';c.fillRect(0,0,W,H);c.font='14px system-ui';
    for(let k=0;k<=4;k++) {
      const y=ymax*k/4;c.strokeStyle='#e4e9ee';c.beginPath();c.moveTo(L,Y(y));c.lineTo(W-R,Y(y));c.stroke();
      c.fillStyle=ink;c.textAlign='right';c.fillText(y.toFixed(ymax>10?0:2),L-8,Y(y)+5);
      const x=xmin+(xmax-xmin)*k/4;c.textAlign='center';c.fillText(x.toFixed(1),X(x),H-B+24);
    }
    c.save();c.beginPath();c.rect(L,T,W-L-R,H-T-B);c.clip();
    for(const s of series) {
      if(s.lo) {
        c.fillStyle=s.fill;c.beginPath();s.x.forEach((x,i)=>i?c.lineTo(X(x),Y(s.lo[i])):c.moveTo(X(x),Y(s.lo[i])));
        for(let i=s.x.length-1;i>=0;i--)c.lineTo(X(s.x[i]),Y(s.hi[i]));c.closePath();c.fill();
      }
      c.strokeStyle=s.color;c.lineWidth=2.5;c.setLineDash(s.dash?[7,5]:[]);c.beginPath();
      s.x.forEach((x,i)=>i?c.lineTo(X(x),Y(s.y[i])):c.moveTo(X(x),Y(s.y[i])));c.stroke();
    }
    c.setLineDash([]);
    for(const p of points){c.beginPath();c.arc(X(p.x),Y(p.y/p.n),4.5,0,2*Math.PI);c.fillStyle=ink;c.fill();}
    c.restore();c.textAlign='center';c.fillStyle=ink;c.fillText(xlabel,(L+W-R)/2,H-5);
    c.save();c.translate(15,H/2);c.rotate(-Math.PI/2);c.fillText(ylabel,0,0);c.restore();
  }
  function density(canvas,x,weights,m,s,name) {
    const step=x[1]-x[0],post=weights.map(w=>w/step);
    const xmin=Math.min(x[0],m-3.5*s),xmax=Math.max(x[x.length-1],m+3.5*s);
    const px=axis(xmin,xmax,501),prior=px.map(v=>normal(v,m,s));
    draw(canvas,[{x:px,y:prior,color:blue,dash:true},{x,y:post,color:red}],
      {xmin,xmax,ymax:Math.max(...post,...prior)*1.12,xlabel:name});
  }
  function setup(root) {
    const logistic=root.dataset.mode==='regression';
    root.innerHTML=`<div class="a09-intro"><strong>${logistic?'Régression logistique : α et β':'Un événement binaire : α → p'}</strong><p>Bleu pointillé : prior · Rouge : postérieur · Intervalles centraux à 90 %.</p></div>
      <div class="a09-controls">
      ${field('ma','Moyenne du prior sur α',logistic?-1:0,-4,4)}
      ${field('sa','Écart-type du prior sur α',1,.25,logistic?3:10,.25)}
      ${logistic?field('mb','Moyenne du prior sur β',0,-3,3)+field('sb','Écart-type du prior sur β',1,.25,3,.25):field('n','Nombre d’essais N',20,0,500,1,'number')+field('y','Nombre de succès y',5,0,500,1,'number')}
      </div>
      ${logistic?`<details class="a09-data"><summary>Changer le jeu de données simulé</summary><p>Ces valeurs servent uniquement à générer les données. Modifier un prior ne régénère pas les observations.</p><div class="a09-controls">${field('truthA','α générateur',-1,-3,3)}${field('truthB','β générateur',1,-3,3)}${field('count','Essais par valeur de x',15,1,100,1,'number')}${field('seed','Graine de simulation',42,1,999999,1,'number')}</div><button type="button" data-action="generate">Générer les données</button><div class="a09-data-table"></div></details>`:''}
      <div class="a09-actions"><button type="button" data-action="reset">Réinitialiser</button><button type="button" data-action="weak">Prior large</button><button type="button" data-action="strong">Prior resserré</button></div>
      <p class="a09-status" role="status" aria-live="polite"></p><div class="a09-results"></div><div class="a09-plots"></div>
      <p class="a09-note">${logistic?'Les bandes décrivent l’incertitude sur p(x), pas un intervalle de futurs résultats binaires. Les points noirs sont les proportions observées.':'N = 0 permet de vérifier que le postérieur est identique au prior. Les probabilités ne sont pas les résultats individuels 0/1.'}</p>
      <details><summary>Méthode de calcul</summary><p>Postérieur proportionnel à la vraisemblance binomiale multipliée par les priors normaux. Intégration numérique sur une grille ${logistic?'bidimensionnelle, puis 2 400 points pondérés pour les bandes de probabilité':'unidimensionnelle'}. Les résultats sont des approximations pédagogiques, sans MCMC. ${logistic?'Les priors sur α et β sont indépendants ; leur postérieur peut être corrélé.':''}</p></details>`;
    const get=k=>Number(root.querySelector(`[data-key="${k}"]`).value);
    const set=(k,v)=>{root.querySelector(`[data-key="${k}"]`).value=v;};
    let data=simulate(-1,1,15,42),timer;
    const status=root.querySelector('.a09-status'),results=root.querySelector('.a09-results'),plots=root.querySelector('.a09-plots');
    const ca=chart(plots,'Distribution de α','Prior et postérieur du coefficient alpha');
    const cb=chart(plots,logistic?'Distribution de β':'Distribution de p','Prior et postérieur');
    const cc=logistic?chart(plots,'Probabilités en fonction de x','Moyenne postérieure et intervalles des probabilités'):null;
    function row(name,st){return `<tr><th scope="row">${name}</th><td>${fmt(st.mean)}</td><td>${fmt(st.sd)}</td><td>${fmt(st.variance)}</td><td>[${fmt(st.lower)}, ${fmt(st.upper)}]</td></tr>`;}
    function showData(){if(logistic)root.querySelector('.a09-data-table').innerHTML='<table><caption>Données actuellement utilisées</caption><tr><th>x</th>'+data.map(d=>`<td>${d.x}</td>`).join('')+'</tr><tr><th>Succès / essais</th>'+data.map(d=>`<td>${d.y}/${d.n}</td>`).join('')+'</tr></table>';}
    function update() {
      for(const input of root.querySelectorAll('input')) input.nextElementSibling.value=input.value;
      if([...root.querySelectorAll('input')].some(input=>!input.checkValidity() || input.value==='')){status.textContent='Veuillez entrer des valeurs dans les limites indiquées.';return;}
      if(!logistic && get('y')>get('n')){status.textContent='Le nombre de succès doit être inférieur ou égal au nombre d’essais.';return;}
      const ma=get('ma'),sa=get('sa');let rows;
      if(logistic) {
        const mb=get('mb'),sb=get('sb'),r=regression(ma,sa,mb,sb,data);
        density(ca,r.a,r.wa,ma,sa,'α');density(cb,r.b,r.wb,mb,sb,'β');
        const x=r.curve.map(d=>d.x);
        draw(cc,[{x,y:r.curve.map(d=>d.priorMedian),lo:r.curve.map(d=>d.priorLo),hi:r.curve.map(d=>d.priorHi),color:blue,fill:'#176caa18',dash:true},
          {x,y:r.curve.map(d=>d.mean),lo:r.curve.map(d=>d.lo),hi:r.curve.map(d=>d.hi),color:red,fill:'#be3d3633'}],{xmin:-2,xmax:2,xlabel:'x standardisé',ylabel:'Probabilité',points:data});
        rows=row('α — postérieur',r.alpha)+row('β — postérieur',r.beta);
        status.textContent=`${sum(data.map(d=>d.y))} succès / ${sum(data.map(d=>d.n))} essais · Corrélation postérieure α–β : ${fmt(r.correlation,2)}.${r.edge>.005?' Attention : masse proche des limites de la grille ; approximation moins fiable.':''}`;
        root._a09Result={...r,data:data.map(d=>({...d}))};
      } else {
        const r=intercept(ma,sa,get('n'),get('y'));density(ca,r.x,r.posterior,ma,sa,'α');
        // Change-of-variable density f_p(p)=f_alpha(logit(p))/(p(1-p)).
        const step=r.x[1]-r.x[0],indices=r.p.map((p,i)=>p>.001&&p<.999?i:-1).filter(i=>i>=0);
        const x=indices.map(i=>r.p[i]);
        const prior=indices.map(i=>r.prior[i]/step/(r.p[i]*(1-r.p[i])));
        const post=indices.map(i=>r.posterior[i]/step/(r.p[i]*(1-r.p[i])));
        draw(cb,[{x,y:prior,color:blue,dash:true},{x,y:post,color:red}],{xmin:0,xmax:1,ymax:Math.max(...prior,...post)*1.1,xlabel:'p = logistique(α)'});
        rows=row('p — prior',r.priorProb)+row('p — postérieur',r.prob)+row('α — postérieur',r.alpha);
        status.textContent=`Données fixes : ${get('y')} succès sur ${get('n')} essais. Densité de p affichée entre 0,001 et 0,999 ; moments calculés sur toute la grille.`;
        root._a09Result=r;
      }
      results.innerHTML='<table><caption>Résumé des distributions</caption><thead><tr><th>Quantité</th><th>Moyenne</th><th>Écart-type</th><th>Variance</th><th>Intervalle à 90 %</th></tr></thead><tbody>'+rows+'</tbody></table>';
      root.dataset.ready='true';
    }
    root.addEventListener('input',event=>{
      if(event.target.tagName!=='INPUT')return;
      event.target.nextElementSibling.value=event.target.value;
      if(['truthA','truthB','count','seed'].includes(event.target.dataset.key))return;
      clearTimeout(timer);timer=setTimeout(update,100);
    });
    root.addEventListener('click',event=>{
      const action=event.target.dataset.action;if(!action)return;
      if(action==='generate') {
        if([...root.querySelectorAll('input')].some(i=>!i.checkValidity() || i.value==='')){status.textContent='Veuillez vérifier les valeurs de simulation.';return;}
        data=simulate(get('truthA'),get('truthB'),get('count'),get('seed'));showData();
      }
      if(action==='weak'){set('sa',logistic?3:10);if(logistic)set('sb',3);}
      if(action==='strong'){set('sa',.5);if(logistic)set('sb',.5);}
      if(action==='reset') {
        for(const input of root.querySelectorAll('input'))input.value=input.defaultValue;
        data=simulate(-1,1,15,42);showData();
      }
      clearTimeout(timer);update();
    });
    showData();update();
  }
  function start(){document.querySelectorAll('.a09-simulator').forEach(setup);}
  if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',start);else start();
})();
