const assert = require('node:assert/strict');
const sim = require('../courses/simulators/a09-priors.js');
const close=(a,b,tol=1e-5)=>assert.ok(Math.abs(a-b)<tol,`${a} != ${b}`);
let r=sim.intercept(0,1,0,0);
close(r.alpha.mean,0);close(r.alpha.variance,1);close(r.prob.mean,.5);
close(r.prob.variance,.043379035858,1e-6);
assert.deepEqual(r.prior,r.posterior);
r=sim.intercept(0,1,20,5);
const mirrored=sim.intercept(0,1,20,15);
close(r.prob.mean+mirrored.prob.mean,1);
assert.ok(r.prob.mean>.25&&r.prob.mean<.5);
assert.ok(sim.intercept(0,1,200,50).prob.variance<r.prob.variance);
assert.ok(sim.intercept(0,.25,20,5).prob.mean>r.prob.mean);
for(const [n,y] of [[0,0],[500,0],[500,500]]) {
 const fit=sim.intercept(-4,10,n,y);
 assert.ok(Number.isFinite(fit.prob.mean)&&fit.prob.mean>0&&fit.prob.mean<1);
}
const data=sim.simulate(-1,1,15,42);
assert.deepEqual(data,sim.simulate(-1,1,15,42));
const empty=sim.regression(-1,1,0,1,[]);
close(empty.alpha.mean,-1);close(empty.alpha.variance,1);close(empty.beta.variance,1);
const post=sim.regression(-1,1,0,1,data);
assert.ok(post.beta.mean>0&&post.beta.sd<1);
const strong=sim.regression(-1,1,0,.25,data);
assert.ok(strong.beta.mean<post.beta.mean);
for(const fit of [post,strong]) for(const p of fit.curve) {
 assert.ok(0<=p.lo&&p.lo<=p.mean&&p.mean<=p.hi&&p.hi<=1);
}
// Check grid convergence including a concentrated likelihood with broad priors.
const big=sim.simulate(-1,1,100,42);
const coarse=sim.regression(3,3,-3,3,big,201);
const fine=sim.regression(3,3,-3,3,big,321);
close(coarse.alpha.mean,fine.alpha.mean,1e-3);
close(coarse.beta.mean,fine.beta.mean,1e-3);
close(coarse.beta.sd,fine.beta.sd,1e-3);
console.log('PASS: prior recovery, symmetry, concentration, regularization, extremes, reproducibility, probability bounds, grid convergence.');
// Section 10: individual Bernoulli data with a reproducible prefix.
const seq20=sim.simulateSequential(-1,1,20,42);
const seq200=sim.simulateSequential(-1,1,200,42);
assert.deepEqual(seq200.slice(0,20),seq20);
assert.notDeepEqual(seq20,sim.simulateSequential(-1,1,20,43));
assert.ok(seq200.every(d=>d.x>=-2&&d.x<=2&&d.n===1&&(d.y===0||d.y===1)));
const early=sim.regression(0,1,0,1,seq20);
const late=sim.regression(0,1,0,1,seq200);
assert.ok(late.alpha.sd<early.alpha.sd&&late.beta.sd<early.beta.sd);
const seqFine=sim.regression(0,1,0,1,seq20,321);
close(early.alpha.mean,seqFine.alpha.mean,1e-3);
close(early.beta.mean,seqFine.beta.mean,1e-3);
// Section 11: use exact odds/probability transformations, including negative changes.
const c=sim.contrast(Math.log(.1/.9),Math.log(2),0,1);
close(c.p0,.1);close(c.p1,2/11);close(c.oddsRatio,2);close(c.difference,2/11-.1);
close(sim.contrast(0,Math.log(2),0,1).oddsRatio,c.oddsRatio);
assert.notEqual(sim.contrast(0,Math.log(2),0,1).difference,c.difference);
close(sim.contrast(-1,0,1,1).difference,0);
close(sim.contrast(-1,1,1,0).difference,0);
assert.ok(sim.contrast(-1,-1,0,1).difference<0);
const reverse=sim.contrast(Math.log(.1/.9),Math.log(2),1,-1);
close(reverse.oddsRatio,.5);close(reverse.difference,-c.difference);
console.log('PASS: sequential prefix/reproducibility, posterior precision, grid convergence, and log-odds/probability interpretation.');
