/* ============================= TOOLTIPS (single shared bubble, portal to body) ============================= */
const sharedTip = document.getElementById('sharedTipBubble');
let pinnedTip = null;

function canHover(){ return window.matchMedia('(hover:hover) and (pointer:fine)').matches; }

function showTip(tip){
  const text = tip.getAttribute('data-tip');
  if(!text) return;
  sharedTip.textContent = text;
  const bw = sharedTip.offsetWidth, bh = sharedTip.offsetHeight; // forces layout with the new text already set
  const rect = tip.getBoundingClientRect();
  let top = rect.top - bh - 10;
  if(top < 8) top = rect.bottom + 10;
  let left = rect.left + rect.width/2 - bw/2;
  left = Math.max(8, Math.min(left, window.innerWidth - bw - 8));
  sharedTip.style.top = top+'px';
  sharedTip.style.left = left+'px';
  sharedTip.classList.add('tip-visible');
}
function hideTip(){ sharedTip.classList.remove('tip-visible'); }

document.addEventListener('mouseover', (e)=>{
  const tip = e.target.closest('.info-tip');
  if(tip && canHover()) showTip(tip);
});
document.addEventListener('mouseout', (e)=>{
  const tip = e.target.closest('.info-tip');
  if(tip && tip!==pinnedTip && (!e.relatedTarget || !tip.contains(e.relatedTarget))) hideTip();
});
document.addEventListener('focusin', (e)=>{
  const tip = e.target.closest('.info-tip');
  if(tip) showTip(tip);
});
document.addEventListener('focusout', (e)=>{
  const tip = e.target.closest('.info-tip');
  if(tip && tip!==pinnedTip) hideTip();
});
document.addEventListener('click', (e)=>{
  const tip = e.target.closest('.info-tip');
  if(tip){
    if(pinnedTip===tip){ pinnedTip=null; hideTip(); }
    else { pinnedTip=tip; showTip(tip); }
  } else if(pinnedTip){ pinnedTip=null; hideTip(); }
});
document.addEventListener('keydown', (e)=>{
  if(e.key==='Escape'){ pinnedTip=null; hideTip(); }
});
window.addEventListener('resize', ()=>{ if(pinnedTip) showTip(pinnedTip); else hideTip(); if(typeof fitFormula==='function') fitFormula(); });
document.addEventListener('scroll', ()=>{ if(!pinnedTip) hideTip(); }, true);

/* ============================= COMPLEX NUMBERS ============================= */
class Cx{
  constructor(re,im=0){this.re=re;this.im=im;}
  add(o){return new Cx(this.re+o.re,this.im+o.im);}
  sub(o){return new Cx(this.re-o.re,this.im-o.im);}
  mul(o){return new Cx(this.re*o.re-this.im*o.im,this.re*o.im+this.im*o.re);}
  div(o){const d=o.re*o.re+o.im*o.im;return new Cx((this.re*o.re+this.im*o.im)/d,(this.im*o.re-this.re*o.im)/d);}
  scale(k){return new Cx(this.re*k,this.im*k);}
  abs(){return Math.hypot(this.re,this.im);}
  arg(){return Math.atan2(this.im,this.re);}
  powInt(n){
    if(n===0) return new Cx(1,0);
    const r=this.abs(), th=this.arg();
    const rn=Math.pow(r,n);
    return new Cx(rn*Math.cos(th*n), rn*Math.sin(th*n));
  }
  isRealish(eps=1e-6){return Math.abs(this.im)<eps;}
}
function clamp(x,a,b){return Math.max(a,Math.min(b,x));}

/* ============================= ROOT SOLVERS ============================= */
function solveQuadratic(c1,c2){
  const b=-c1, c=-c2;
  const disc=b*b-4*c;
  if(disc>=0){
    const sq=Math.sqrt(disc);
    return [new Cx((-b+sq)/2),new Cx((-b-sq)/2)];
  } else {
    const sq=Math.sqrt(-disc);
    return [new Cx(-b/2,sq/2),new Cx(-b/2,-sq/2)];
  }
}
function solveCubic(c1,c2,c3){
  const p=-c1, q=-c2, r=-c3;
  const P=q-(p*p)/3;
  const Q=(2*p*p*p)/27-(p*q)/3+r;
  const shift=-p/3;
  const disc=(Q*Q)/4+(P*P*P)/27;
  let ys=[];
  if(Math.abs(P)<1e-9 && Math.abs(Q)<1e-9){
    ys=[new Cx(0),new Cx(0),new Cx(0)];
  } else if(Math.abs(disc)<1e-9){
    const u=Math.cbrt(-Q/2);
    ys=[new Cx(2*u),new Cx(-u),new Cx(-u)];
  } else if(disc>0){
    const sq=Math.sqrt(disc);
    const u=Math.cbrt(-Q/2+sq);
    const v=Math.cbrt(-Q/2-sq);
    const y1=u+v;
    const reP=-(u+v)/2;
    const imP=(u-v)*Math.sqrt(3)/2;
    ys=[new Cx(y1),new Cx(reP,imP),new Cx(reP,-imP)];
  } else {
    const m=2*Math.sqrt(-P/3);
    const theta=Math.acos(clamp((3*Q)/(P*m),-1,1))/3;
    for(let k=0;k<3;k++){
      ys.push(new Cx(m*Math.cos(theta-(2*Math.PI*k)/3)));
    }
  }
  return ys.map(y=>y.add(new Cx(shift)));
}


function groupRoots(roots,tol=1e-4){
  const groups=[];
  for(const rt of roots){
    const found=groups.find(g=>g.root.sub(rt).abs()<tol*(1+g.root.abs()));
    if(found) found.mult++; else groups.push({root:rt,mult:1});
  }
  return groups;
}

/* ============================= LINEAR SOLVER (complex) ============================= */
function solveComplexSystem(M,b){
  const n=b.length;
  const A=M.map(row=>row.map(v=>v));
  const rhs=b.map(v=>v);
  for(let col=0;col<n;col++){
    let maxRow=col,maxVal=A[col][col].abs();
    for(let r=col+1;r<n;r++){ if(A[r][col].abs()>maxVal){maxVal=A[r][col].abs();maxRow=r;} }
    if(maxRow!==col){ [A[col],A[maxRow]]=[A[maxRow],A[col]]; [rhs[col],rhs[maxRow]]=[rhs[maxRow],rhs[col]]; }
    const pivot=A[col][col];
    for(let r=col+1;r<n;r++){
      const factor=A[r][col].div(pivot);
      for(let c=col;c<n;c++){ A[r][c]=A[r][c].sub(factor.mul(A[col][c])); }
      rhs[r]=rhs[r].sub(factor.mul(rhs[col]));
    }
  }
  const x=new Array(n);
  for(let i=n-1;i>=0;i--){
    let sum=rhs[i];
    for(let j=i+1;j<n;j++){ sum=sum.sub(A[i][j].mul(x[j])); }
    x[i]=sum.div(A[i][i]);
  }
  return x;
}

/* ============================= EXPLICIT MODEL BUILDERS ============================= */
function buildHomogeneousModel(order, initials, coeffs){
  const roots = order===2 ? solveQuadratic(coeffs[0],coeffs[1]) : solveCubic(coeffs[0],coeffs[1],coeffs[2]);
  const groups = groupRoots(roots);
  const basisMeta=[];
  groups.forEach((g,gi)=>{ for(let j=0;j<g.mult;j++){ basisMeta.push({gi,j,root:g.root}); } });
  const k=basisMeta.length;
  const M=[];
  for(let n=0;n<k;n++){
    const row=[];
    for(const bm of basisMeta){ row.push(bm.root.powInt(n).scale(Math.pow(n,bm.j))); }
    M.push(row);
  }
  const b = initials.map(v=>new Cx(v,0));
  const coeffsSolved = solveComplexSystem(M,b);
  return {groups, basisMeta, coeffs:coeffsSolved};
}

function evalHomogeneous(model,n){
  let total=new Cx(0);
  model.basisMeta.forEach((bm,i)=>{
    total = total.add(model.coeffs[i].mul(bm.root.powInt(n)).scale(Math.pow(n,bm.j)));
  });
  return total.re;
}
function order1Explicit(n,a0,c1,c0){
  if(Math.abs(c1-1)<1e-9) return a0 + c0*n;
  return Math.pow(c1,n)*a0 + c0*(Math.pow(c1,n)-1)/(c1-1);
}


function fmtNum(x,d=4){
  if(Object.is(x,-0)) x=0;
  const r=Math.round(x*Math.pow(10,d))/Math.pow(10,d);
  return r.toString();
}


function toFraction(x, maxDen=24, eps=1e-6){
  if(Math.abs(x) < eps) return {num:0, den:1};
  const sign = x<0 ? -1 : 1; x = Math.abs(x);
  let h1=1,h2=0,k1=0,k2=1,b=x,a,aux,iter=0;
  do{
    a = Math.floor(b);
    aux=h1; h1=a*h1+h2; h2=aux;
    aux=k1; k1=a*k1+k2; k2=aux;
    if(Math.abs(b-a) < 1e-9) break;
    b = 1/(b-a);
    iter++;
  } while(k1<=maxDen && iter<20);
  if(k1>maxDen || k1===0) return null;
  if(Math.abs(x - h1/k1) > eps*Math.max(1,x)) return null;
  return {num:sign*h1, den:k1};
}
function latexNum(x){
  if(Math.abs(x - Math.round(x)) < 1e-7) return `${Math.round(x)}`;
  const f = toFraction(x);
  if(f && f.den>1){
    const sign = f.num<0 ? '-' : '';
    return `${sign}\\frac{${Math.abs(f.num)}}{${f.den}}`;
  }
  return fmtNum(x,4);
}
function latexPow(baseLatex,n){ return `\\left(${baseLatex}\\right)^{${n}}`; }

function buildLatexHomogeneous(model){
  const parts=[]; const seen=new Set();
  model.groups.forEach((g,gi)=>{
    if(seen.has(gi)) return;
    const coeffIdxStart = model.basisMeta.findIndex(bm=>bm.gi===gi && bm.j===0);
    if(g.root.isRealish()){
      const rL = latexNum(g.root.re); let poly;
      if(g.mult===1){ poly = latexNum(model.coeffs[coeffIdxStart].re); }
      else if(g.mult===2){
        const bIdx=model.basisMeta.findIndex(bm=>bm.gi===gi && bm.j===1);
        poly = `\\left(${latexNum(model.coeffs[coeffIdxStart].re)} + ${latexNum(model.coeffs[bIdx].re)}n\\right)`;
      } else {
        const bIdx=model.basisMeta.findIndex(bm=>bm.gi===gi && bm.j===1);
        const cIdx=model.basisMeta.findIndex(bm=>bm.gi===gi && bm.j===2);
        poly = `\\left(${latexNum(model.coeffs[coeffIdxStart].re)} + ${latexNum(model.coeffs[bIdx].re)}n + ${latexNum(model.coeffs[cIdx].re)}n^2\\right)`;
      }
      parts.push(`${poly}${latexPow(rL,'n')}`);
      seen.add(gi);
    } else {
      const partnerIdx = model.groups.findIndex((h,hi)=> !seen.has(hi) && hi!==gi && Math.abs(h.root.re-g.root.re)<1e-4 && Math.abs(h.root.im+g.root.im)<1e-4);
      const c = model.coeffs[coeffIdxStart];
      const r = g.root.abs(); const theta = g.root.arg();
      const alpha = 2*c.re, beta = -2*c.im;
      parts.push(`${latexPow(latexNum(r),'n')}\\left(${latexNum(alpha)}\\cos(${latexNum(theta)}n) + ${latexNum(beta)}\\sin(${latexNum(theta)}n)\\right)`);
      seen.add(gi); if(partnerIdx>=0) seen.add(partnerIdx);
    }
  });
  const body = parts.join(' + ').replace(/\+ -/g,'- ').replace(/\+ \\left\(-/g,'- \\left(');
  return `a_n = ${body} \\quad ; \\quad n \\geq 0`;
}
function buildLatexOrder1(a0,c1,c0){
  if(Math.abs(c1-1)<1e-9) return `a_n = ${latexNum(a0)} + ${latexNum(c0)}n \\quad ; \\quad n \\geq 0`;
  const fixed = c0/(1-c1);
  const sign = fixed<0 ? '-' : '+';
  return `a_n = \\left(${latexNum(a0-fixed)}\\right)${latexPow(latexNum(c1),'n')} ${sign} ${latexNum(Math.abs(fixed))} \\quad ; \\quad n \\geq 0`;
}

/* ============================= PROCESS / "VER PROCESO" HELPERS ============================= */
function recurrenceLatex(order, coeffs){
  let str='';
  coeffs.forEach((c,i)=>{
    const term = `${latexNum(Math.abs(c))}\\,a_{n-${i+1}}`;
    str += (i===0 ? (c<0?'-':'') : (c<0?' - ':' + ')) + term;
  });
  return `a_n = ${str}`;
}
function charEqLatex(order, coeffs){
  const terms=[`t^{${order}}`];
  for(let i=0;i<order;i++){
    const c=coeffs[i], power=order-1-i;
    const sign = c>=0 ? '-' : '+';
    const absVal = Math.abs(c);
    const coefStr = (Math.abs(absVal-1)<1e-9 && power>0) ? '' : latexNum(absVal);
    const varPart = power===0 ? '' : (power===1 ? 't' : `t^{${power}}`);
    terms.push(`${sign} ${coefStr}${varPart}`);
  }
  return terms.join(' ') + ' = 0';
}
function describeRootNature(model){
  const allReal = model.groups.every(g=>g.root.isRealish());
  if(allReal){
    if(model.groups.length===model.basisMeta.length) return 'reales y distintas';
    if(model.groups.length===1) return `real repetida (multiplicidad ${model.groups[0].mult})`;
    return 'reales, con al menos una repetida';
  }
  const hasReal = model.groups.some(g=>g.root.isRealish());
  return hasReal ? 'una raíz real y un par de raíces complejas conjugadas' : 'complejas conjugadas';
}
function rootsListLatex(model){
  return model.groups.map(g=>{
    if(g.root.isRealish()) return `t = ${latexNum(g.root.re)}` + (g.mult>1 ? ` \\text{ (multiplicidad ${g.mult})}` : '');
    return `t = ${latexNum(g.root.re)} ${g.root.im>=0?'+':'-'} ${latexNum(Math.abs(g.root.im))}i`;
  }).join(' \\quad ; \\quad ');
}
function buildGeneralFormLatex(model){
  const pool=['\\alpha','\\beta','\\gamma','\\delta','\\epsilon'];
  let li=0; const next=()=> pool[li++] || `k_{${li}}`;
  const parts=[]; const seen=new Set();
  model.groups.forEach((g,gi)=>{
    if(seen.has(gi)) return;
    if(g.root.isRealish()){
      const rL=latexNum(g.root.re); let poly;
      if(g.mult===1){ poly=next(); }
      else if(g.mult===2){ const a=next(),b=next(); poly=`\\left(${a} + ${b}\\,n\\right)`; }
      else { const a=next(),b=next(),c=next(); poly=`\\left(${a} + ${b}\\,n + ${c}\\,n^2\\right)`; }
      parts.push(`${poly}${latexPow(rL,'n')}`);
      seen.add(gi);
    } else {
      const partnerIdx = model.groups.findIndex((h,hi)=> !seen.has(hi) && hi!==gi && Math.abs(h.root.re-g.root.re)<1e-4 && Math.abs(h.root.im+g.root.im)<1e-4);
      const rL=latexNum(g.root.abs());
      const a=next(), b=next();
      parts.push(`${latexPow(rL,'n')}\\left(${a}\\cos(n\\theta) + ${b}\\sin(n\\theta)\\right)`);
      seen.add(gi); if(partnerIdx>=0) seen.add(partnerIdx);
    }
  });
  return `a_n = ${parts.join(' + ')} \\quad ; \\quad n \\geq 0`;
}
function buildInitialConditionEquations(model, initials){
  const pool=['\\alpha','\\beta','\\gamma'];
  const k = model.basisMeta.length;
  const rows=[];
  for(let n=0;n<k;n++){
    const termStrs=[];
    model.basisMeta.forEach((bm,i)=>{
      const rVal = bm.root.re;
      const coefNumeric = Math.pow(rVal,n)*Math.pow(n,bm.j);
      if(Math.abs(coefNumeric) < 1e-12) return;
      let coefStr;
      if(Math.abs(coefNumeric-1)<1e-9) coefStr='';
      else if(Math.abs(coefNumeric+1)<1e-9) coefStr='-';
      else coefStr = latexNum(coefNumeric);
      termStrs.push(`${coefStr}${pool[i]}`);
    });
    const rhs = termStrs.join(' + ').replace(/\+ -/g,'- ');
    rows.push(`a_{${n}} = ${rhs} = ${latexNum(initials[n])}`);
  }
  return rows;
}
function buildSolvedConstants(model){
  const pool=['\\alpha','\\beta','\\gamma'];
  return model.basisMeta.map((bm,i)=> `${pool[i]} = ${latexNum(model.coeffs[i].re)}`).join(' \\quad ; \\quad ');
}
function buildOrder1Steps(a0,c1,c0){
  const steps=[];
  steps.push({title:'1. Partimos de la relación recursiva', eq:`a_n = ${latexNum(c1)}\\,a_{n-1} + ${latexNum(c0)} \\quad ; \\quad a_0 = ${latexNum(a0)}`});
  if(Math.abs(c1-1)<1e-9){
    steps.push({title:'2. El multiplicador es 1', eq:`a_n = a_{n-1} + ${latexNum(c0)}`, note:'Cada semana solo se suma una cantidad fija: es una progresión aritmética.'});
    steps.push({title:'3. Fórmula explícita', eq: buildLatexOrder1(a0,c1,c0)});
  } else {
    const fixed = c0/(1-c1);
    steps.push({title:'2. Buscamos el punto fijo', eq:`x^{*} = \\dfrac{${latexNum(c0)}}{1-${latexNum(c1)}} = ${latexNum(fixed)}`, note:'Es el valor donde aₙ dejaría de cambiar (aₙ = aₙ₋₁).'});
    steps.push({title:'3. Restamos el punto fijo a ambos lados', eq:`a_n - x^{*} = ${latexNum(c1)}\\left(a_{n-1} - x^{*}\\right)`, note:`Esto convierte la recurrencia en una progresión geométrica pura, de razón ${latexNum(c1)}.`});
    steps.push({title:'4. Resolviendo la geométrica', eq:`a_n - x^{*} = (${latexNum(c1)})^n\\left(a_0 - x^{*}\\right)`});
    steps.push({title:'5. Fórmula explícita', eq: buildLatexOrder1(a0,c1,c0)});
  }
  return steps;
}
function buildComprobacion(order, coeffs, recSeq, evalExplicitFn, k){
  const rv = recSeq[k];
  let recEqLatex, recNote='';
  if(k < order){
    recEqLatex = `a_{${k}} = ${latexNum(rv)}`;
    recNote = 'Condición inicial dada (no se calcula con la recurrencia).';
  } else {
    const terms=[], valTerms=[];
    for(let j=1;j<=order;j++){
      terms.push(`${latexNum(coeffs[j-1])}\\,a_{${k-j}}`);
      valTerms.push(`${latexNum(coeffs[j-1])}(${latexNum(recSeq[k-j])})`);
    }
    recEqLatex = `a_{${k}} = ${terms.join(' + ').replace(/\+ -/g,'- ')} = ${valTerms.join(' + ').replace(/\+ -/g,'- ')} = ${latexNum(rv)}`;
  }
  const ev = evalExplicitFn(k);
  return {recEqLatex, recNote, expEqLatex:`a_{${k}} = ${latexNum(ev)}`, match: Math.abs(rv-ev)<1e-4};
}

/* ============================= SEQUENCE GENERATION ============================= */
function generateIterative(order, initials, coeffs, n){
  const arr=new Array(n);
  if(order===1){
    arr[0]=initials[0];
    for(let i=1;i<n;i++) arr[i]=coeffs[0]*arr[i-1]+coeffs[1];
  } else if(order===2){
    arr[0]=initials[0]; arr[1]=initials[1];
    for(let i=2;i<n;i++) arr[i]=coeffs[0]*arr[i-1]+coeffs[1]*arr[i-2];
  } else {
    arr[0]=initials[0]; arr[1]=initials[1]; arr[2]=initials[2];
    for(let i=3;i<n;i++) arr[i]=coeffs[0]*arr[i-1]+coeffs[1]*arr[i-2]+coeffs[2]*arr[i-3];
  }
  return arr;
}
function naiveRec(order, initials, coeffs, n){
  if(order===1){
    const f=(m)=> m===0 ? initials[0] : coeffs[0]*f(m-1)+coeffs[1];
    return f(n);
  } else if(order===2){
    const f=(m)=> m===0?initials[0] : m===1?initials[1] : coeffs[0]*f(m-1)+coeffs[1]*f(m-2);
    return f(n);
  } else {
    const f=(m)=> m===0?initials[0] : m===1?initials[1] : m===2?initials[2] : coeffs[0]*f(m-1)+coeffs[1]*f(m-2)+coeffs[2]*f(m-3);
    return f(n);
  }
}
function estimateNaiveCalls(order, benchN){
  if(order===1) return benchN+1;
  const base = order===2 ? 1.7 : 2.2;
  return Math.pow(base, benchN);
}

/* ============================= PLAIN-LANGUAGE STATUS ============================= */
function classifyTrend(seq){
  const n = seq.length;
  const diffs = [];
  for(let i=1;i<n;i++) diffs.push(seq[i]-seq[i-1]);
  let signChanges=0;
  for(let i=1;i<diffs.length;i++){
    const a=Math.sign(diffs[i-1]), b=Math.sign(diffs[i]);
    if(a!==0 && b!==0 && a!==b) signChanges++;
  }
  const oscillating = diffs.length>=4 && signChanges >= Math.max(2, Math.floor(diffs.length*0.35));

  const first=seq[0], last=seq[n-1];
  const base = Math.abs(first)<1e-9 ? 1 : Math.abs(first);
  const pct = ((last-first)/base)*100;

  if(oscillating){
    return {css:'status-osc', icon:'🔁', label:'Más o menos — demanda oscilante',
      msg:`Sube y baja de forma repetida a lo largo de las ${n} semanas. No hay una sola tendencia: conviene planear para picos y valles, no para un promedio.`};
  }
  if(pct > 25) return {css:'status-up', icon:'📈', label:'Todo va bien — crecimiento fuerte',
    msg:`La demanda proyectada crece ${fmtNum(pct,1)}% entre la semana 0 y la ${n-1}. Va a hacer falta más capacidad pronto.`};
  if(pct > 3) return {css:'status-up', icon:'📈', label:'Todo va bien — crecimiento sano',
    msg:`La demanda proyectada crece ${fmtNum(pct,1)}% en ${n} semanas, a un ritmo manejable.`};
  if(pct >= -3) return {css:'status-flat', icon:'➖', label:'Más o menos — demanda estable',
    msg:`Se mantiene casi igual (${fmtNum(pct,1)}%) a lo largo de las ${n} semanas.`};
  if(pct >= -25) return {css:'status-down', icon:'📉', label:'Cuidado — demanda a la baja',
    msg:`La demanda proyectada cae ${fmtNum(Math.abs(pct),1)}% en ${n} semanas.`};
  return {css:'status-down', icon:'📉', label:'Todo va mal — caída fuerte',
    msg:`La demanda proyectada cae ${fmtNum(Math.abs(pct),1)}% en ${n} semanas. Conviene revisar el servicio.`};
}

/* ============================= UI STATE ============================= */
let currentOrder=1;
let currentRunContext=null;
let chartInstance=null;

const defaults = {
  1:{initials:[980], coeffs:[1.04,15], labels:['a₀ (semana base)'], coefLabels:['c₁ (multiplicador)','c₀ (carga fija)'],
     tips:['Dato ya observado, no se calcula: la demanda real que tuviste esa semana.'],
     coefTips:['Multiplica la demanda de la semana anterior (ej. 1.04 = crece 4% cada semana).','Cantidad fija que se suma cada semana, además del crecimiento.']},
  2:{initials:[1200,1260], coeffs:[1.05,-0.02], labels:['a₀','a₁'], coefLabels:['c₁','c₂'],
     tips:['Demanda observada en la semana 0.','Demanda observada en la semana 1.'],
     coefTips:['Peso de la demanda de hace 1 semana sobre la actual.','Peso de la demanda de hace 2 semanas sobre la actual.']},
  3:{initials:[800,840,905], coeffs:[0.6,0.3,0.08], labels:['a₀','a₁','a₂'], coefLabels:['c₁','c₂','c₃'],
     tips:['Demanda observada en la semana 0.','Demanda observada en la semana 1.','Demanda observada en la semana 2.'],
     coefTips:['Peso de hace 1 semana.','Peso de hace 2 semanas.','Peso de hace 3 semanas.']}
};

function tipSpan(text){
  return `<span class="info-tip" tabindex="0" role="button" aria-label="Ayuda" data-tip="${text}">ⓘ</span>`;
}

function renderDynamicFields(order){
  const d = defaults[order];
  const container = document.getElementById('dynamicFields');
  let html = '<div class="field-row">';
  d.labels.forEach((lab,i)=>{
    html += `<div class="field-group field-enter" style="animation-delay:${i*40}ms;"><label for="init${i}">${lab}${tipSpan(d.tips[i])}</label><input type="number" step="any" id="init${i}" value="${d.initials[i]}"></div>`;
  });
  html += '</div><div class="field-row">';
  d.coefLabels.forEach((lab,i)=>{
    const delay=(d.labels.length+i)*40;
    html += `<div class="field-group field-enter" style="animation-delay:${delay}ms;"><label for="coef${i}">${lab}${tipSpan(d.coefTips[i])}</label><input type="number" step="any" id="coef${i}" value="${d.coeffs[i]}"></div>`;
  });
  html += '</div>';
  container.innerHTML = html;

  const benchInput = document.getElementById('benchTerm');
  if(order===1) benchInput.max=60;
  else if(order===2) benchInput.max=30;
  else benchInput.max=20;
  if(parseInt(benchInput.value,10) > parseInt(benchInput.max,10)) benchInput.value = benchInput.max;
}

document.querySelectorAll('.order-toggle button').forEach((btn,idx)=>{
  btn.addEventListener('click', ()=>{
    document.querySelectorAll('.order-toggle button').forEach(b=>b.setAttribute('aria-pressed','false'));
    btn.setAttribute('aria-pressed','true');
    document.querySelector('.order-toggle').style.setProperty('--idx', idx);
    currentOrder = parseInt(btn.dataset.order,10);
    renderDynamicFields(currentOrder);
  });
});
renderDynamicFields(1);

function readInputs(order){
  const d = defaults[order];
  const initials = d.labels.map((_,i)=>parseFloat(document.getElementById('init'+i).value));
  const coeffs = d.coefLabels.map((_,i)=>parseFloat(document.getElementById('coef'+i).value));
  return {initials, coeffs};
}
function showWarn(msg){
  const w = document.getElementById('warnMsg');
  if(!msg){ w.style.display='none'; w.textContent=''; return; }
  w.style.display='block'; w.textContent = msg;
}
function replay(el){ el.classList.remove('reveal-block'); void el.offsetWidth; el.classList.add('reveal-block'); }
function replayBadge(el){ el.style.animation='none'; void el.offsetWidth; el.style.animation=''; }

function fitFormula(){
  const box = document.getElementById('formulaOut');
  box.style.setProperty('--katex-scale', 1);
  const avail = box.clientWidth;
  const natural = box.scrollWidth;
  if(natural > avail && avail > 0){
    box.style.setProperty('--katex-scale', Math.max(0.55, avail/natural));
  }
}

document.getElementById('runBtn').addEventListener('click', run);

let isolated = false;
document.getElementById('isolateBtn').addEventListener('click', ()=>{
  if(!chartInstance) return;
  isolated = !isolated;
  chartInstance.setDatasetVisibility(0, !isolated);
  chartInstance.update();
  document.getElementById('isolateBtn').textContent = isolated ? 'Mostrar ambas' : 'Aislar explícita';
});

function run(){
  showWarn('');
  isolated = false;
  document.getElementById('isolateBtn').textContent = 'Aislar explícita';
  const order = currentOrder;
  const {initials, coeffs} = readInputs(order);
  let n = parseInt(document.getElementById('nInput').value,10);
  n = clamp(n, order+2, 120);
  document.getElementById('nInput').value = n;

  if(order===2 && Math.abs(coeffs[1])<1e-9){ showWarn('c₂ no puede ser 0 en una recurrencia de orden 2.'); return; }
  if(order===3 && Math.abs(coeffs[2])<1e-9){ showWarn('c₃ no puede ser 0 en una recurrencia de orden 3.'); return; }

  const recSeq = generateIterative(order, initials, coeffs, n);

  let evalExplicit, formulaLatex;
  if(order===1){
    evalExplicit = (m)=> order1Explicit(m, initials[0], coeffs[0], coeffs[1]);
    formulaLatex = buildLatexOrder1(initials[0], coeffs[0], coeffs[1]);
  } else {
    const model = buildHomogeneousModel(order, initials, coeffs);
    evalExplicit = (m)=> evalHomogeneous(model, m);
    formulaLatex = buildLatexHomogeneous(model);
  }
  katex.render(formulaLatex, document.getElementById('formulaOut'), {throwOnError:false, displayMode:true});
  fitFormula();

  const status = classifyTrend(recSeq);
  const statusCard = document.getElementById('statusCard');
  statusCard.className = 'status-card reveal-block ' + status.css;
  document.getElementById('statusIcon').textContent = status.icon;
  document.getElementById('statusLabel').textContent = status.label;
  document.getElementById('statusMsg').textContent = status.msg;

  const tbody = document.getElementById('termsBody');
  tbody.innerHTML='';
  let maxDiff=0; const expSeq=[];
  for(let i=0;i<n;i++){
    const rv = recSeq[i], ev = evalExplicit(i);
    expSeq.push(ev);
    const diff = Math.abs(rv-ev);
    maxDiff = Math.max(maxDiff, diff);
    const tr = document.createElement('tr');
    tr.innerHTML = `<td>${i}<button type="button" class="proc-btn" data-n="${i}">Ver Proceso</button></td><td>${fmtNum(rv,4)}</td><td>${fmtNum(ev,4)}</td><td class="${diff<1e-6?'diff-ok':'diff-bad'}">${diff<1e-9?'0':fmtNum(diff,6)}</td>`;
    tbody.appendChild(tr);
  }

  currentRunContext = {
    order, initials, coeffs,
    model: (order===1 ? null : buildHomogeneousModel(order, initials, coeffs)),
    recSeq, formulaLatex
  };
  const badge = document.getElementById('verifyBadge');
  const detail = document.getElementById('verifyDetail');
  badge.style.display='inline-flex';
  if(maxDiff<1e-4){ badge.className='badge ok'; badge.textContent=`✓ Coinciden en los ${n} términos`; }
  else { badge.className='badge warn'; badge.textContent='⚠ Divergencia detectada'; }
  detail.textContent = `diferencia máxima: ${fmtNum(maxDiff,6)}`;
  replayBadge(badge);

  const ctx = document.getElementById('seqChart').getContext('2d');
  if(chartInstance) chartInstance.destroy();
  chartInstance = new Chart(ctx,{
    type:'line',
    data:{
      labels:Array.from({length:n},(_,i)=>i),
      datasets:[
        {label:'Recurrente', data:recSeq, borderColor:'#e3a34b', backgroundColor:'rgba(227,163,75,0.12)', pointRadius:2, borderWidth:2.5, tension:0.15},
        {label:'Explícita', data:expSeq, borderColor:'#2f8f86', backgroundColor:'rgba(47,143,134,0.1)', pointRadius:0, borderWidth:2, borderDash:[5,3], tension:0.15}
      ]
    },
    options:{
      responsive:true,
      animation:{duration:500,easing:'easeOutQuint'},
      plugins:{legend:{position:'top', labels:{boxWidth:14, color:'#171526'}}},
      scales:{ x:{title:{display:true,text:'Semana (n)'}}, y:{title:{display:true,text:'Demanda proyectada'}} }
    }
  });

  // ---- Benchmark with runaway-recursion safety guard ----
  let benchN = clamp(parseInt(document.getElementById('benchTerm').value,10), order, parseInt(document.getElementById('benchTerm').max,10));
  let reps = clamp(parseInt(document.getElementById('repsInput').value,10), 1, 5000);

  const PER_CALL_CEILING = 2000000;
  while(estimateNaiveCalls(order,benchN) > PER_CALL_CEILING && benchN>order){ benchN--; }
  const estCalls = estimateNaiveCalls(order,benchN);
  const SAFE_TOTAL = 3000000;
  let adjusted = false;
  if(estCalls*reps > SAFE_TOTAL){
    const newReps = Math.max(1, Math.floor(SAFE_TOTAL/estCalls));
    if(newReps < reps){ reps = newReps; adjusted = true; }
  }
  document.getElementById('benchTerm').value = benchN;
  document.getElementById('repsInput').value = reps;
  if(adjusted){ showWarn(`Se ajustaron término/repeticiones para no bloquear el navegador (la recursión pura crece exponencialmente aquí).`); }

  const t0 = performance.now();
  let sinkA=0;
  for(let i=0;i<reps;i++){ sinkA += naiveRec(order, initials, coeffs, benchN); }
  const t1 = performance.now();
  let sinkB=0;
  for(let i=0;i<reps;i++){ sinkB += evalExplicit(benchN); }
  const t2 = performance.now();

  const recTime = t1-t0, expTime = t2-t1;
  document.getElementById('benchRec').textContent = recTime.toFixed(2)+' ms';
  document.getElementById('benchExp').textContent = expTime.toFixed(3)+' ms';
  const maxTime = Math.max(recTime, expTime, 0.001);
  requestAnimationFrame(()=>{
    document.getElementById('barRec').style.transform = `scaleX(${recTime/maxTime})`;
    document.getElementById('barExp').style.transform = `scaleX(${expTime/maxTime})`;
  });

  const speedup = expTime>0 ? (recTime/expTime) : (recTime>0 ? Infinity : 1);
  document.getElementById('speedupNote').textContent = isFinite(speedup)
    ? `Con ${reps} repeticiones sobre a(${benchN}): la fórmula explícita fue ${speedup>=1 ? speedup.toFixed(1)+'× más rápida' : (1/speedup).toFixed(1)+'× más lenta'} que la recursión pura.`
    : `La fórmula explícita resolvió a(${benchN}) en un tiempo no medible frente a la recursión pura.`;

  ['formulaOut','chartWrap'].length; // no-op guard
  [document.querySelector('.formula-box'), ...document.querySelectorAll('.reveal-block')].forEach(el=>{ if(el) replay(el); });
}

/* ============================= PROCESS MODAL ============================= */
function addModalStep(body, title, eqLatex, note){
  const h = document.createElement('h4'); h.textContent = title; body.appendChild(h);
  if(eqLatex){
    const div = document.createElement('div'); div.className='step-eq';
    katex.render(eqLatex, div, {throwOnError:false, displayMode:true});
    body.appendChild(div);
  }
  if(note){ const p=document.createElement('p'); p.className='step-note'; p.textContent=note; body.appendChild(p); }
}

function openProcessModal(k){
  if(!currentRunContext) return;
  const {order, initials, coeffs, model, recSeq} = currentRunContext;
  const evalExplicitFn = order===1
    ? (m)=> order1Explicit(m, initials[0], coeffs[0], coeffs[1])
    : (m)=> evalHomogeneous(model, m);

  document.getElementById('modalTitle').textContent = `Proceso — término a\u2099 para n = ${k}`;
  const body = document.getElementById('modalBody');
  body.innerHTML = '';

  if(order===1){
    buildOrder1Steps(initials[0], coeffs[0], coeffs[1]).forEach(s=> addModalStep(body, s.title, s.eq, s.note));
  } else {
    addModalStep(body, '1. Sustituimos aₙ = tⁿ en la relación', recurrenceLatex(order,coeffs).replace('a_n =', 't^{n} ='));
    addModalStep(body, '2. Hacemos cero un lado y factorizamos', charEqLatex(order,coeffs), 'Esta es la ecuación característica asociada a la sucesión.');
    const nature = describeRootNature(model);
    addModalStep(body, '3. Resolvemos la ecuación característica', rootsListLatex(model), `Las raíces son ${nature}.`);
    addModalStep(body, '4. Forma general de la solución', buildGeneralFormLatex(model));
    const allReal = model.groups.every(g=>g.root.isRealish());
    if(allReal){
      const rows = buildInitialConditionEquations(model, initials);
      rows.forEach((r,i)=> addModalStep(body, i===0 ? '5. Usamos las condiciones iniciales' : '', r));
      addModalStep(body, '6. Resolviendo el sistema', buildSolvedConstants(model));
      addModalStep(body, '7. Fórmula explícita', currentRunContext.formulaLatex);
    } else {
      addModalStep(body, '5. Ajustamos las constantes', buildSolvedConstants(model), 'Se resuelven planteando el sistema con las condiciones iniciales (a₀, a₁, ...) e igualando parte real e imaginaria.');
      addModalStep(body, '6. Fórmula explícita', currentRunContext.formulaLatex);
    }
  }

  const comp = buildComprobacion(order, coeffs, recSeq, evalExplicitFn, k);
  const compTitle = document.createElement('h4');
  compTitle.textContent = `Comprobación para n = ${k}`;
  body.appendChild(compTitle);
  const grid = document.createElement('div'); grid.className='compare-grid';
  const col1 = document.createElement('div'); col1.className='compare-col';
  col1.innerHTML = '<div class="cc-label">Recurrente</div>';
  const eq1 = document.createElement('div'); katex.render(comp.recEqLatex, eq1, {throwOnError:false});
  col1.appendChild(eq1);
  if(comp.recNote){ const p=document.createElement('p'); p.className='step-note'; p.textContent=comp.recNote; col1.appendChild(p); }
  const col2 = document.createElement('div'); col2.className='compare-col';
  col2.innerHTML = '<div class="cc-label">Explícita</div>';
  const eq2 = document.createElement('div'); katex.render(comp.expEqLatex, eq2, {throwOnError:false});
  col2.appendChild(eq2);
  grid.appendChild(col1); grid.appendChild(col2);
  body.appendChild(grid);
  const verdict = document.createElement('p'); verdict.className='step-note';
  verdict.style.marginTop='10px'; verdict.style.fontWeight='600';
  verdict.style.color = comp.match ? 'var(--turq)' : 'var(--danger)';
  verdict.textContent = comp.match ? '✓ Ambos métodos coinciden.' : '⚠ Los métodos no coinciden — revisa los datos ingresados.';
  body.appendChild(verdict);

  lastFocusedEl = document.activeElement;
  const overlay = document.getElementById('modalOverlay');
  overlay.classList.add('open');
  document.getElementById('modalClose').focus();
  document.addEventListener('keydown', onModalKeydown);
}
function closeProcessModal(){
  document.getElementById('modalOverlay').classList.remove('open');
  document.removeEventListener('keydown', onModalKeydown);
  if(lastFocusedEl && typeof lastFocusedEl.focus==='function') lastFocusedEl.focus();
}
function onModalKeydown(e){ if(e.key==='Escape') closeProcessModal(); }
let lastFocusedEl=null;

document.getElementById('termsBody').addEventListener('click', (e)=>{
  const btn = e.target.closest('.proc-btn');
  if(btn) openProcessModal(parseInt(btn.dataset.n,10));
});
document.getElementById('modalClose').addEventListener('click', closeProcessModal);
document.getElementById('modalOverlay').addEventListener('click', (e)=>{
  if(e.target.id==='modalOverlay') closeProcessModal();
});

run();
