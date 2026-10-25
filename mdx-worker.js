/* نواة فصل الصوت (MDX-Net) — بدون اعتماديات. تعمل داخل Worker أو Node */
(function (G) {
  'use strict';
  /* ---------- FFT مختلط الأساس (2/3/4/5) على نمط kissfft ---------- */
  function FFT(n) {
    this.n = n;
    const tw = (this.tw = new Float64Array(2 * n));
    for (let i = 0; i < n; i++) { const a = -2 * Math.PI * i / n; tw[2 * i] = Math.cos(a); tw[2 * i + 1] = Math.sin(a); }
    const f = []; let m = n;
    for (const p of [4, 2, 3, 5]) while (m % p === 0) { f.push(p, m / p); m /= p; }
    if (m !== 1) throw new Error('FFT size unsupported: ' + n);
    this.f = f;
    this.sr = new Float64Array(n); this.si = new Float64Array(n);
  }
  FFT.prototype._work = function (oR, oI, o, iR, iI, io, fs, is, fi) {
    const p = this.f[fi], m = this.f[fi + 1], tw = this.tw, end = o + p * m;
    if (m === 1) { for (let q = o; q < end; q++) { oR[q] = iR[io]; oI[q] = iI[io]; io += fs * is; } }
    else { for (let q = o; q < end; q += m) { this._work(oR, oI, q, iR, iI, io, fs * p, is, fi + 2); io += fs * is; } }
    const n = this.n;
    if (p === 2) {
      for (let k = 0; k < m; k++) {
        const a = o + k, b = a + m, t = 2 * k * fs, wr = tw[t], wi = tw[t + 1];
        const tr = oR[b] * wr - oI[b] * wi, ti = oR[b] * wi + oI[b] * wr;
        oR[b] = oR[a] - tr; oI[b] = oI[a] - ti; oR[a] += tr; oI[a] += ti;
      }
    } else if (p === 4) {
      for (let k = 0; k < m; k++) {
        const a = o + k, b = a + m, c = b + m, d = c + m, t = 2 * k * fs;
        const w1r = tw[t], w1i = tw[t + 1], w2r = tw[2 * t], w2i = tw[2 * t + 1], w3r = tw[3 * t], w3i = tw[3 * t + 1];
        const s0r = oR[b] * w1r - oI[b] * w1i, s0i = oR[b] * w1i + oI[b] * w1r;
        const s1r = oR[c] * w2r - oI[c] * w2i, s1i = oR[c] * w2i + oI[c] * w2r;
        const s2r = oR[d] * w3r - oI[d] * w3i, s2i = oR[d] * w3i + oI[d] * w3r;
        const s5r = oR[a] - s1r, s5i = oI[a] - s1i;
        oR[a] += s1r; oI[a] += s1i;
        const s3r = s0r + s2r, s3i = s0i + s2i, s4r = s0r - s2r, s4i = s0i - s2i;
        oR[c] = oR[a] - s3r; oI[c] = oI[a] - s3i; oR[a] += s3r; oI[a] += s3i;
        oR[b] = s5r + s4i; oI[b] = s5i - s4r; oR[d] = s5r - s4i; oI[d] = s5i + s4r;
      }
    } else if (p === 3) {
      const ei = tw[2 * fs * m + 1];
      for (let k = 0; k < m; k++) {
        const a = o + k, b = a + m, c = b + m, t = 2 * k * fs;
        const w1r = tw[t], w1i = tw[t + 1], w2r = tw[2 * t], w2i = tw[2 * t + 1];
        const s1r = oR[b] * w1r - oI[b] * w1i, s1i = oR[b] * w1i + oI[b] * w1r;
        const s2r = oR[c] * w2r - oI[c] * w2i, s2i = oR[c] * w2i + oI[c] * w2r;
        const s3r = s1r + s2r, s3i = s1i + s2i; let s0r = (s1r - s2r) * ei, s0i = (s1i - s2i) * ei;
        oR[b] = oR[a] - s3r * 0.5; oI[b] = oI[a] - s3i * 0.5;
        oR[a] += s3r; oI[a] += s3i;
        oR[c] = oR[b] + s0i; oI[c] = oI[b] - s0r; oR[b] -= s0i; oI[b] += s0r;
      }
    } else if (p === 5) {
      const yar = tw[2 * fs * m], yai = tw[2 * fs * m + 1], ybr = tw[4 * fs * m], ybi = tw[4 * fs * m + 1];
      for (let u = 0; u < m; u++) {
        const i0 = o + u, i1 = i0 + m, i2 = i1 + m, i3 = i2 + m, i4 = i3 + m, t = 2 * u * fs;
        const s0r = oR[i0], s0i = oI[i0];
        let wr = tw[t], wi = tw[t + 1]; const s1r = oR[i1] * wr - oI[i1] * wi, s1i = oR[i1] * wi + oI[i1] * wr;
        wr = tw[2 * t]; wi = tw[2 * t + 1]; const s2r = oR[i2] * wr - oI[i2] * wi, s2i = oR[i2] * wi + oI[i2] * wr;
        wr = tw[3 * t]; wi = tw[3 * t + 1]; const s3r = oR[i3] * wr - oI[i3] * wi, s3i = oR[i3] * wi + oI[i3] * wr;
        wr = tw[4 * t]; wi = tw[4 * t + 1]; const s4r = oR[i4] * wr - oI[i4] * wi, s4i = oR[i4] * wi + oI[i4] * wr;
        const s7r = s1r + s4r, s7i = s1i + s4i, s10r = s1r - s4r, s10i = s1i - s4i;
        const s8r = s2r + s3r, s8i = s2i + s3i, s9r = s2r - s3r, s9i = s2i - s3i;
        oR[i0] = s0r + s7r + s8r; oI[i0] = s0i + s7i + s8i;
        const s5r = s0r + s7r * yar + s8r * ybr, s5i = s0i + s7i * yar + s8i * ybr;
        const s6r = s10i * yai + s9i * ybi, s6i = -s10r * yai - s9r * ybi;
        oR[i1] = s5r - s6r; oI[i1] = s5i - s6i; oR[i4] = s5r + s6r; oI[i4] = s5i + s6i;
        const s11r = s0r + s7r * ybr + s8r * yar, s11i = s0i + s7i * ybr + s8i * yar;
        const s12r = -s10i * ybi + s9i * yai, s12i = s10r * ybi - s9r * yai;
        oR[i2] = s11r + s12r; oI[i2] = s11i + s12i; oR[i3] = s11r - s12r; oI[i3] = s11i - s12i;
      }
    }
  };
  /* تحويل مباشر: (iR,iI) → (oR,oI) */
  FFT.prototype.fwd = function (iR, iI, oR, oI) { this._work(oR, oI, 0, iR, iI, 0, 1, 1, 0); };
  /* عكسي (مع القسمة على n) عبر تبديل الجزأين */
  FFT.prototype.inv = function (iR, iI, oR, oI) {
    this._work(oI, oR, 0, iI, iR, 0, 1, 1, 0);
    const s = 1 / this.n; for (let i = 0; i < this.n; i++) { oR[i] *= s; oI[i] *= s; }
  };

  /* ---------- أدوات STFT لقناتين (L في الحقيقي، R في التخيلي) ---------- */
  function Spec(P) {
    const n = P.nfft; this.P = P; this.fft = new FFT(n);
    this.win = new Float64Array(n); for (let i = 0; i < n; i++) this.win[i] = 0.5 - 0.5 * Math.cos(2 * Math.PI * i / n);
    this.zr = new Float64Array(n); this.zi = new Float64Array(n); this.yr = new Float64Array(n); this.yi = new Float64Array(n);
    const hop = P.hop, T = P.dimT, len = n + hop * (T - 1);
    this.env = new Float64Array(len);
    for (let t = 0; t < T; t++) for (let i = 0; i < n; i++) this.env[t * hop + i] += this.win[i] * this.win[i];
    this.chunk = hop * (T - 1);
  }
  const refl = (i, N) => { if (i < 0) i = -i; if (i >= N) i = 2 * (N - 1) - i; return i < 0 ? 0 : i; };
  /* L,R: Float32Array بطول chunk → Float32Array [4*dimF*dimT] */
  Spec.prototype.stft = function (L, R, out) {
    const P = this.P, n = P.nfft, hop = P.hop, T = P.dimT, F = P.dimF, N = L.length, h = n >> 1, zr = this.zr, zi = this.zi, yr = this.yr, yi = this.yi, w = this.win;
    for (let t = 0; t < T; t++) {
      const s = t * hop - h;
      for (let i = 0; i < n; i++) { const j = refl(s + i, N); zr[i] = L[j] * w[i]; zi[i] = R[j] * w[i]; }
      this.fft.fwd(zr, zi, yr, yi);
      for (let k = 0; k < F; k++) {
        const kk = k === 0 ? 0 : n - k;
        const lr = (yr[k] + yr[kk]) * .5, li = (yi[k] - yi[kk]) * .5, rr = (yi[k] + yi[kk]) * .5, ri = -(yr[k] - yr[kk]) * .5;
        out[(0 * F + k) * T + t] = lr; out[(1 * F + k) * T + t] = li; out[(2 * F + k) * T + t] = rr; out[(3 * F + k) * T + t] = ri;
      }
    }
  };
  /* spec [4*dimF*dimT] → outL,outR بطول chunk */
  Spec.prototype.istft = function (spec, outL, outR) {
    const P = this.P, n = P.nfft, hop = P.hop, T = P.dimT, F = P.dimF, h = n >> 1, zr = this.zr, zi = this.zi, yr = this.yr, yi = this.yi, w = this.win, env = this.env;
    const len = env.length, accL = new Float64Array(len), accR = new Float64Array(len);
    for (let t = 0; t < T; t++) {
      zr.fill(0); zi.fill(0);
      for (let k = 0; k < F && k <= h; k++) {
        let ar = spec[(0 * F + k) * T + t], ai = spec[(1 * F + k) * T + t], br = spec[(2 * F + k) * T + t], bi = spec[(3 * F + k) * T + t];
        if (k === 0 || k === h) { ai = 0; bi = 0; }
        zr[k] = ar - bi; zi[k] = ai + br;
        if (k !== 0 && k !== h) { zr[n - k] = ar + bi; zi[n - k] = -ai + br; }
      }
      this.fft.inv(zr, zi, yr, yi);
      const s = t * hop;
      for (let i = 0; i < n; i++) { accL[s + i] += yr[i] * w[i]; accR[s + i] += yi[i] * w[i]; }
    }
    const C = this.chunk;
    for (let i = 0; i < C; i++) { const e = env[i + h]; const g = e > 1e-11 ? 1 / e : 0; outL[i] = accL[i + h] * g; outR[i] = accR[i + h] * g; }
  };

  /* ---------- الفصل الكامل لمسار واحد ---------- */
  /* P: {nfft,dimF,dimT,hop,compensate}; run(Float32Array input)->Promise<Float32Array> */
  async function separate(P, run, L, R, onProg, shouldStop) {
    const sp = new Spec(P), chunk = sp.chunk, trim = P.nfft >> 1, gen = chunk - 2 * trim, N = L.length;
    const pad = gen - (N % gen || gen) % gen;
    const total = trim + N + pad + trim;
    const mL = new Float32Array(total), mR = new Float32Array(total);
    mL.set(L, trim); mR.set(R, trim);
    const oL = new Float32Array(N), oR = new Float32Array(N);
    const inp = new Float32Array(4 * P.dimF * P.dimT), cl = new Float32Array(chunk), cr = new Float32Array(chunk), tl = new Float32Array(chunk), tr = new Float32Array(chunk);
    const steps = Math.ceil((N + pad) / gen);
    for (let s = 0, i = 0; i < N + pad; i += gen, s++) {
      if (shouldStop && shouldStop()) throw new Error('__cx');
      cl.set(mL.subarray(i, i + chunk)); cr.set(mR.subarray(i, i + chunk));
      sp.stft(cl, cr, inp);
      for (let c = 0; c < 4; c++) for (let k = 0; k < 3; k++) inp.fill(0, (c * P.dimF + k) * P.dimT, (c * P.dimF + k + 1) * P.dimT);
      const out = await run(inp);
      sp.istft(out, tl, tr);
      const m = Math.min(gen, N - i), cmp = P.compensate || 1;
      for (let j = 0; j < m; j++) { oL[i + j] = tl[trim + j] * cmp; oR[i + j] = tr[trim + j] * cmp; }
      if (onProg) onProg((s + 1) / steps);
    }
    return [oL, oR];
  }
  G.MDX = { FFT, Spec, separate };
  if (typeof module !== 'undefined') module.exports = G.MDX;
})(typeof self !== 'undefined' ? self : globalThis);
/* استخراج الأصوات المحيطة/الطبيعية (ضوضاء الخلفية الثابتة: رياح، مطر، مياه، زحام) من المسار المتبقي.
   الفكرة: أرضية الضوضاء لكل تردد = نسبة دنيا زمنية داخل نافذة ~2 ثانية، ثم تنعيمها عبر التردد (وسيط) لتتجاهل القمم النغمية للآلات. */
(function (G) {
  'use strict';
  const N = 2048, H = 512, K = N / 2 + 1, BF = 172, P = (G.AMBP || {}), PCT = P.pct || 0.1, GAIN = P.gain || 3.0, SPW = P.spw || 24, SPQ = P.spq == null ? 0.15 : P.spq;
  function ambient(L, R, onp) {
    const fft = new G.MDX.FFT(N), len = L.length, h = N >> 1;
    const win = new Float64Array(N); for (let i = 0; i < N; i++) win[i] = 0.5 - 0.5 * Math.cos(2 * Math.PI * i / N);
    const T = Math.ceil((len + N) / H) + 1, zr = new Float64Array(N), zi = new Float64Array(N), yr = new Float64Array(N), yi = new Float64Array(N);
    const frame = t => { const s = t * H - h; for (let i = 0; i < N; i++) { const j = s + i; if (j >= 0 && j < len) { zr[i] = L[j] * win[i]; zi[i] = R[j] * win[i]; } else { zr[i] = 0; zi[i] = 0; } } fft.fwd(zr, zi, yr, yi); };
    const mid = k => { const kk = k === 0 ? 0 : N - k; const lr = (yr[k] + yr[kk]) * .5, li = (yi[k] - yi[kk]) * .5, rr = (yi[k] + yi[kk]) * .5, ri = -(yr[k] - yr[kk]) * .5; return Math.sqrt((lr * lr + li * li + rr * rr + ri * ri) * .5); };
    /* مرحلة 1: أرضية لكل كتلة */
    const NB = Math.max(1, Math.ceil(T / BF)), floors = new Float32Array(NB * K), blk = new Float32Array(BF * K), col = new Float32Array(BF), sm = new Float32Array(K), win2 = new Float32Array(2 * SPW + 1);
    for (let b = 0; b < NB; b++) {
      const t0 = b * BF, t1 = Math.min(T, t0 + BF), nf = t1 - t0;
      for (let t = t0; t < t1; t++) { frame(t); for (let k = 0; k < K; k++) blk[(t - t0) * K + k] = mid(k); }
      const q = Math.min(nf - 1, Math.floor(PCT * nf));
      for (let k = 0; k < K; k++) { for (let t = 0; t < nf; t++) col[t] = blk[t * K + k]; const c = col.subarray(0, nf).sort(); sm[k] = c[q]; }
      for (let k = 0; k < K; k++) { let c = 0; for (let j = -SPW; j <= SPW; j++) { const kk = Math.min(K - 1, Math.max(0, k + j)); win2[c++] = sm[kk]; } floors[b * K + k] = win2.slice().sort()[Math.floor(SPQ * 2 * SPW)]; }
      if (onp) onp(.5 * (b + 1) / NB);
    }
    /* مرحلة 2: قناع ناعم + إعادة بناء */
    const OL = new Float64Array(len + N * 2), OR = new Float64Array(len + N * 2), env = new Float64Array(len + N * 2), fl = new Float32Array(K);
    for (let t = 0; t < T; t++) {
      frame(t);
      const bc = (t + .5) / BF - .5, b0 = Math.max(0, Math.min(NB - 1, Math.floor(bc))), b1 = Math.min(NB - 1, b0 + 1), w1 = Math.max(0, Math.min(1, bc - b0));
      for (let k = 0; k < K; k++) fl[k] = floors[b0 * K + k] * (1 - w1) + floors[b1 * K + k] * w1;
      const ar = zr, ai = zi; /* نعيد استخدام المصفوفات لبناء الطيف المقنّع */
      const Yr = Float64Array.from(yr), Yi = Float64Array.from(yi);
      zr.fill(0); zi.fill(0);
      for (let k = 0; k < K; k++) {
        const kk = k === 0 ? 0 : N - k;
        let lr = (Yr[k] + Yr[kk]) * .5, li = (Yi[k] - Yi[kk]) * .5, rr = (Yi[k] + Yi[kk]) * .5, ri = -(Yr[k] - Yr[kk]) * .5;
        const m = Math.sqrt((lr * lr + li * li + rr * rr + ri * ri) * .5), g = m > 1e-9 ? Math.min(1, GAIN * fl[k] / m) : 1;
        lr *= g; li *= g; rr *= g; ri *= g; if (k === 0 || k === h) { li = 0; ri = 0; }
        zr[k] = lr - ri; zi[k] = li + rr;
        if (k !== 0 && k !== h) { zr[N - k] = lr + ri; zi[N - k] = -li + rr; }
      }
      fft.inv(zr, zi, yr, yi);
      const s = t * H - h;
      for (let i = 0; i < N; i++) { const j = s + i; if (j >= 0 && j < len) { OL[j] += yr[i] * win[i]; OR[j] += yi[i] * win[i]; env[j] += win[i] * win[i]; } }
      if (onp && (t & 255) === 0) onp(.5 + .5 * t / T);
    }
    const aL = new Float32Array(len), aR = new Float32Array(len);
    for (let i = 0; i < len; i++) { const e = env[i] > 1e-6 ? 1 / env[i] : 0; aL[i] = OL[i] * e; aR[i] = OR[i] * e; }
    return [aL, aR];
  }
  G.AMB = { ambient };
  if (typeof module !== 'undefined') module.exports = G.AMB;
})(typeof self !== 'undefined' ? self : globalThis);
/* ---- Worker: يعتمد على ort (onnxruntime-web) و MDX المحمّلين قبله ---- */
(function(){
  const S=self;let sess=null,mix=null,stopF=false,job=0,paused=false,pst=0,pms=0;
  const post=(m,t)=>S.postMessage(m,t||[]);
  const tick=async()=>{while(paused)await new Promise(r=>setTimeout(r,120));await new Promise(r=>setTimeout(r,0))};
  S.onmessage=async e=>{
    const m=e.data;
    try{
      if(m.t==='init'){
        ort.env.wasm.wasmPaths=m.paths||m.base;ort.env.wasm.simd=true;ort.env.wasm.proxy=false;
        ort.env.wasm.numThreads=(S.crossOriginIsolated&&m.threads>1)?m.threads:1;
        post({t:'init',iso:!!S.crossOriginIsolated,threads:ort.env.wasm.numThreads});
      }else if(m.t==='load'){
        if(sess){try{await sess.release()}catch(_){}sess=null}
        const bytes=new Uint8Array(m.buf);
        try{sess=await ort.InferenceSession.create(bytes,{executionProviders:['wasm'],graphOptimizationLevel:'all'})}
        catch(er){ if(ort.env.wasm.numThreads>1){ort.env.wasm.numThreads=1;sess=await ort.InferenceSession.create(bytes,{executionProviders:['wasm'],graphOptimizationLevel:'all'})}else throw er}
        post({t:'loaded',id:m.id,threads:ort.env.wasm.numThreads});
      }else if(m.t==='mix'){
        mix=[new Float32Array(m.L),new Float32Array(m.R)];post({t:'mix'});
      }else if(m.t==='sep'){
        stopF=false;pms=0;const P=m.P;let n=0,t0=Date.now();
        const run=async inp=>{await tick();const r=await sess.run({input:new ort.Tensor('float32',inp,[1,4,P.dimF,P.dimT])});n++;return r.output.data};
        const res=await MDX.separate(P,run,mix[0],mix[1],p=>post({t:'p',job:m.job,p,el:(Date.now()-t0-pms)/1000}),()=>stopF);
        post({t:'done',job:m.job,L:res[0],R:res[1],el:(Date.now()-t0-pms)/1000},[res[0].buffer,res[1].buffer]);
      }else if(m.t==='amb'){const r=AMB.ambient(new Float32Array(m.L),new Float32Array(m.R),p=>post({t:'p',job:m.job,p,el:0}));post({t:'ambdone',L:r[0],R:r[1]},[r[0].buffer,r[1].buffer])}
      else if(m.t==='pause'){if(m.v){if(!paused){paused=true;pst=Date.now()}}else if(paused){pms+=Date.now()-pst;paused=false}}
      else if(m.t==='stop'){stopF=true}
      else if(m.t==='free'){mix=null;if(sess){try{await sess.release()}catch(_){}sess=null}post({t:'freed'})}
    }catch(er){post({t:'err',job:m.job,e:String((er&&er.message)||er)})}
  };
  post({t:'ready'});
})();
