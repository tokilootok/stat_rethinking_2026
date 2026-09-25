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
  function simulateSequential(a,b,n,seed) {
    let state=seed>>>0;
    const random=()=>{state=(Math.imul(1664525,state)+1013904223)>>>0;return state/4294967296;};
    return Array.from({length:n},()=>{
      const x=-2+4*random();
      return {x,n:1,y:random()<sigmoid(a+b*x)?1:0};
    });
  }
  function contrast(alpha,beta,x0,dx) {
    const eta0=alpha+beta*x0,eta1=alpha+beta*(x0+dx);
    const p0=sigmoid(eta0),p1=sigmoid(eta1);
    return {eta0,eta1,p0,p1,logOddsChange:beta*dx,oddsRatio:Math.exp(beta*dx),difference:p1-p0,riskRatio:p1/p0};
  }
  const api={sigmoid,normal,stats,intercept,regression,simulate,simulateSequential,contrast};
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
  function draw(canvas,series,{xmin,xmax,ymin=0,ymax=1,xlabel='',ylabel='Densité',points=[],references=[]}) {
    const c=canvas.getContext('2d'),W=760,H=290,L=62,R=20,T=12,B=49;
    const X=x=>L+(x-xmin)/(xmax-xmin)*(W-L-R),Y=y=>H-B-(y-ymin)/(ymax-ymin)*(H-T-B);
    c.clearRect(0,0,W,H);c.fillStyle='#fff';c.fillRect(0,0,W,H);c.font='14px system-ui';
    for(let k=0;k<=4;k++) {
      const y=ymin+(ymax-ymin)*k/4;c.strokeStyle='#e4e9ee';c.beginPath();c.moveTo(L,Y(y));c.lineTo(W-R,Y(y));c.stroke();
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
    c.setLineDash([4,3]);
    for(const x of references){c.strokeStyle='#15804a';c.lineWidth=2;c.beginPath();c.moveTo(X(x),T);c.lineTo(X(x),H-B);c.stroke();}
    c.setLineDash([]);
    for(const p of points){c.beginPath();c.arc(X(p.x),Y(p.y/p.n),4.5,0,2*Math.PI);c.fillStyle=ink;c.fill();}
    c.restore();c.textAlign='center';c.fillStyle=ink;c.fillText(xlabel,(L+W-R)/2,H-5);
    c.save();c.translate(15,H/2);c.rotate(-Math.PI/2);c.fillText(ylabel,0,0);c.restore();
  }
  function density(canvas,x,weights,m,s,name,truth=null) {
    const step=x[1]-x[0],post=weights.map(w=>w/step);
    const xmin=Math.min(x[0],m-3.5*s),xmax=Math.max(x[x.length-1],m+3.5*s);
    const px=axis(xmin,xmax,501),prior=px.map(v=>normal(v,m,s));
    draw(canvas,[{x:px,y:prior,color:blue,dash:true},{x,y:post,color:red}],
      {xmin,xmax,ymax:Math.max(...post,...prior)*1.12,xlabel:name,references:truth===null?[]:[truth]});
  }
  function setupSequential(root) {
    root.innerHTML=`<div class="a09-intro"><strong>Apprendre une courbe logistique, observation après observation</strong><p>Vert : vérité simulée · Bleu pointillé : prior · Rouge : postérieur · Bandes à 90 %.</p></div>
      <fieldset class="a09-group"><legend>1. Choisir le monde à simuler</legend><div class="a09-controls">
        ${field('truthA','α vrai (générateur)',-1,-4,4)}${field('truthB','β vrai (générateur)',1,-4,4)}
        ${field('count','Nombre total de points',20,1,200,1,'number')}${field('seed','Graine de simulation',42,1,999999,1,'number')}
      </div><p>Les x sont tirés uniformément entre −2 et 2. Chaque point est un résultat binaire 0 ou 1.</p>
      <button type="button" data-action="generate">Générer les données</button><p class="a09-pending" role="status"></p></fieldset>
      <fieldset class="a09-group"><legend>2. Choisir les croyances de départ</legend><div class="a09-controls">
        ${field('ma','Moyenne du prior sur α',0,-4,4)}${field('sa','Écart-type du prior sur α',1,.25,3,.25)}
        ${field('mb','Moyenne du prior sur β',0,-4,4)}${field('sb','Écart-type du prior sur β',1,.25,3,.25)}
      </div><p>Modifier les priors recalcule le postérieur sur les mêmes observations déjà révélées. Les valeurs vraies ne sont pas transmises au modèle.</p></fieldset>
      <fieldset class="a09-group"><legend>3. Révéler les observations</legend>
        ${field('seen','Nombre de points révélés',0,0,20,1)}
        <div class="a09-actions"><button type="button" data-action="play">Lancer l’animation</button><button type="button" data-action="next">+1 point</button><button type="button" data-action="all">Tout révéler</button><button type="button" data-action="rewind">Revenir au prior</button><button type="button" data-action="reset">Réinitialiser</button></div>
      </fieldset>
      <p class="a09-status" role="status" aria-live="polite"></p><div class="a09-results"></div><div class="a09-plots"></div>
      <p>À zéro point, le rouge et le bleu se superposent. Les bandes portent sur la probabilité p(x), pas sur les futurs résultats individuels 0/1. La vérité est un repère de simulation, pas une contrainte imposée au postérieur.</p>
      <details><summary>Comment le postérieur est-il mis à jour ?</summary><p>À l’étape t, la grille intègre le prior initial multiplié par la vraisemblance des t premiers points, chacun compté une seule fois. Revenir en arrière retire les derniers points ; changer les priors conserve les données. L’animation est une approximation numérique sur grille, et non un calcul MCMC ou quap. Les bandes sont pointwise (à chaque x), pas une bande simultanée sur toute la courbe.</p></details>`;
    const el=k=>root.querySelector(`[data-key="${k}"]`),get=k=>Number(el(k).value);
    const status=root.querySelector('.a09-status'),pending=root.querySelector('.a09-pending');
    const results=root.querySelector('.a09-results'),plots=root.querySelector('.a09-plots');
    const cp=chart(plots,'La courbe vraie et ce que les données permettent d’apprendre','Courbe vraie verte, prior bleu, postérieur rouge, observations binaires noires');
    const ca=chart(plots,'α : position de la courbe','Prior et postérieur de alpha avec repère de la valeur vraie');
    const cb=chart(plots,'β : pente de la courbe','Prior et postérieur de beta avec repère de la valeur vraie');
    let active={a:-1,b:1,n:20,seed:42},data=simulateSequential(-1,1,20,42),seen=0,playing=false,timer=null,debounce=null;
    const play=root.querySelector('[data-action="play"]');
    function pause(){playing=false;clearTimeout(timer);play.textContent='Lancer l’animation';play.setAttribute('aria-pressed','false');}
    function valid(keys){return keys.every(k=>el(k).value!==''&&el(k).checkValidity());}
    function sync(){for(const input of root.querySelectorAll('input'))input.nextElementSibling.value=input.value;}
    function dirty(){return get('truthA')!==active.a||get('truthB')!==active.b||get('count')!==active.n||get('seed')!==active.seed;}
    function showPending(){pending.textContent=dirty()?'Paramètres de génération modifiés : cliquez sur « Générer les données » pour les appliquer. Les graphes montrent encore la simulation précédente.':'';}
    function update(){
      sync();
      if(!valid(['ma','sa','mb','sb','seen'])){pause();status.textContent='Veuillez vérifier les valeurs des priors et du nombre de points révélés.';return false;}
      seen=get('seen');const observed=data.slice(0,seen);
      const r=regression(get('ma'),get('sa'),get('mb'),get('sb'),observed);
      const x=r.curve.map(d=>d.x);
      draw(cp,[
        {x,y:r.curve.map(d=>d.priorMedian),lo:r.curve.map(d=>d.priorLo),hi:r.curve.map(d=>d.priorHi),color:blue,fill:'#176caa18',dash:true},
        {x,y:r.curve.map(d=>d.mean),lo:r.curve.map(d=>d.lo),hi:r.curve.map(d=>d.hi),color:red,fill:'#be3d3633'},
        {x,y:x.map(v=>sigmoid(active.a+active.b*v)),color:'#15804a'}
      ],{xmin:-2,xmax:2,xlabel:'x (entre −2 et 2)',ylabel:'Probabilité / résultat',points:observed});
      density(ca,r.a,r.wa,get('ma'),get('sa'),'α',active.a);
      density(cb,r.b,r.wb,get('mb'),get('sb'),'β',active.b);
      const row=(name,trueValue,st)=>`<tr><th scope="row">${name}</th><td>${fmt(trueValue,2)}</td><td>${fmt(st.mean)}</td><td>${fmt(st.sd)}</td><td>[${fmt(st.lower)}, ${fmt(st.upper)}]</td></tr>`;
      results.innerHTML='<table><caption>Paramètres vrais et estimations à cette étape</caption><thead><tr><th>Paramètre</th><th>Vrai</th><th>Moyenne postérieure</th><th>Écart-type</th><th>Intervalle à 90 %</th></tr></thead><tbody>'+row('α',active.a,r.alpha)+row('β',active.b,r.beta)+'</tbody></table>';
      const last=seen?` Dernier point : x = ${fmt(observed[seen-1].x,2)}, y = ${observed[seen-1].y}.`:'';
      status.textContent=`${seen} / ${data.length} points révélés · ${sum(observed.map(d=>d.y))} succès · Corrélation α–β : ${fmt(r.correlation,2)}.${last}${r.edge>.005?' Masse proche du bord de grille : approximation moins fiable.':''}`;
      root._a09Result={...r,data:observed.map(d=>({...d})),allData:data.map(d=>({...d})),seen,active:{...active}};
      root.dataset.ready='true';root.dataset.seen=String(seen);
      root.querySelector('[data-action="next"]').disabled=seen>=data.length;
      showPending();return true;
    }
    function tick(){
      if(!playing)return;
      el('seen').value=Math.min(data.length,get('seen')+1);
      if(!update()||get('seen')>=data.length){pause();return;}
      timer=setTimeout(tick,650);
    }
    function generate(){
      pause();clearTimeout(debounce);
      if(!valid(['truthA','truthB','count','seed'])){pending.textContent='Veuillez entrer des paramètres valides (1 à 200 points, graine entière positive).';return;}
      active={a:get('truthA'),b:get('truthB'),n:get('count'),seed:get('seed')};
      data=simulateSequential(active.a,active.b,active.n,active.seed);
      el('seen').max=active.n;el('seen').value=0;update();
    }
    root.addEventListener('input',event=>{
      if(event.target.tagName!=='INPUT')return;
      pause();clearTimeout(debounce);event.target.nextElementSibling.value=event.target.value;
      if(['truthA','truthB','count','seed'].includes(event.target.dataset.key)){showPending();return;}
      debounce=setTimeout(update,100);
    });
    root.addEventListener('click',event=>{
      const action=event.target.dataset.action;if(!action)return;
      clearTimeout(debounce);
      if(action==='play'){
        if(playing){pause();return;}
        if(get('seen')>=data.length)el('seen').value=0;
        if(!update())return;
        playing=true;play.textContent='Pause';play.setAttribute('aria-pressed','true');timer=setTimeout(tick,100);return;
      }
      pause();
      if(action==='generate'){generate();return;}
      if(action==='reset'){for(const input of root.querySelectorAll('input'))input.value=input.defaultValue;generate();return;}
      if(action==='next')el('seen').value=Math.min(data.length,get('seen')+1);
      if(action==='all')el('seen').value=data.length;
      if(action==='rewind')el('seen').value=0;
      update();
    });
    document.addEventListener('visibilitychange',()=>{if(document.hidden)pause();});
    update();
  }

  function setupInterpretation(root) {
    root.innerHTML=`<div class="a09-intro"><strong>Interpréter β : qu’est-ce qui change vraiment ?</strong><p>Ici α et β sont des valeurs choisies pour illustrer le modèle, pas des priors ni des estimations postérieures.</p></div>
      <div class="a09-controls">${field('alpha','α : intercept',-2.1972245773362196,-5,5,'any')}${field('beta','β : coefficient',0.6931471805599453,-4,4,'any')}${field('x0','x₀ : point de départ',0,-2,2,.1)}${field('dx','Δx : changement du prédicteur',1,-2,2,.1)}</div>
      <div class="a09-actions"><button type="button" data-action="double">Exemple : odds ×2, départ à 10 %</button><button type="button" data-action="null">β = 0</button><button type="button" data-action="negative">Inverser le signe de β</button></div>
      <p class="a09-status" role="status" aria-live="polite"></p><div class="a09-results"></div><div class="a09-plots"></div>
      <p>La courbe est calculée avec p(x) = logistique(α + βx). Les points noirs repèrent les probabilités des scénarios comparés ; ce ne sont pas des observations 0/1. Aucun effet causal n’est établi par ce calcul.</p>`;
    const get=k=>Number(root.querySelector(`[data-key="${k}"]`).value);
    const set=(k,v)=>{root.querySelector(`[data-key="${k}"]`).value=v;};
    const plots=root.querySelector('.a09-plots');
    const curve=chart(plots,'Deux scénarios sur la même courbe logistique','Probabilités aux points de départ et d’arrivée');
    const difference=chart(plots,'Le même odds ratio à différents niveaux de départ','Différence de probabilité selon la probabilité initiale');
    function update(){
      for(const input of root.querySelectorAll('input'))input.nextElementSibling.value=fmt(input.value,3);
      const a=get('alpha'),b=get('beta'),x0=get('x0'),dx=get('dx'),r=contrast(a,b,x0,dx);
      const pc=v=>fmt(100*v,2)+' %';
      root.querySelector('.a09-status').textContent=`De x₀ = ${fmt(x0,2)} à x₁ = ${fmt(x0+dx,2)} : βΔx = ${fmt(r.logOddsChange)} ; odds multipliées par ${fmt(r.oddsRatio)}.`;
      root.querySelector('.a09-results').innerHTML=`<table><caption>Suivre le calcul d’une échelle à l’autre</caption><thead><tr><th>Échelle</th><th>Départ x₀</th><th>Arrivée x₁</th><th>Comparaison</th></tr></thead><tbody>
        <tr><th>Log-odds : α + βx</th><td>${fmt(r.eta0)}</td><td>${fmt(r.eta1)}</td><td>Différence : ${fmt(r.logOddsChange)}</td></tr>
        <tr><th>Odds : exp(α + βx)</th><td>${fmt(Math.exp(r.eta0))}</td><td>${fmt(Math.exp(r.eta1))}</td><td>Rapport : ${fmt(r.oddsRatio)}</td></tr>
        <tr><th>Probabilité : logistique(α + βx)</th><td>${pc(r.p0)}</td><td>${pc(r.p1)}</td><td>${fmt(100*r.difference,2)} points de pourcentage</td></tr>
        <tr><th>Rapport de probabilités</th><td colspan="2">p₁ / p₀</td><td>${fmt(r.riskRatio)} (variation relative : ${fmt(100*(r.riskRatio-1),2)} %)</td></tr></tbody></table>`;
      const x=axis(-4.04,4.04,101);
      draw(curve,[{x,y:x.map(v=>sigmoid(a+b*v)),color:blue}],{xmin:-4,xmax:4,xlabel:'x',ylabel:'Probabilité',points:[{x:x0,y:r.p0,n:1},{x:x0+dx,y:r.p1,n:1}]});
      const p=axis(0,1,151),delta=p.map(v=>r.oddsRatio*v/(1-v+r.oddsRatio*v)-v);
      draw(difference,[{x:p,y:delta,color:red}],{xmin:0,xmax:1,ymin:Math.min(0,...delta)*1.12,ymax:Math.max(.01,...delta)*1.12,xlabel:'Probabilité de départ p₀',ylabel:'Différence p₁ − p₀',points:[{x:r.p0,y:r.difference,n:1}]});
      root._a09Result=r;root.dataset.ready='true';
    }
    root.addEventListener('input',update);
    root.addEventListener('click',event=>{
      const action=event.target.dataset.action;if(!action)return;
      if(action==='double'){set('alpha',Math.log(.1/.9));set('beta',Math.log(2));set('x0',0);set('dx',1);}
      if(action==='null')set('beta',0);
      if(action==='negative')set('beta',-get('beta'));
      update();
    });
    update();
  }

  function setup(root) {
    if(root.dataset.mode==='interpretation'){setupInterpretation(root);return;}
    if(root.dataset.mode==='sequential'){setupSequential(root);return;}
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
