// Cloudflare Worker — وسيط يسمح لنسخة الويب بالاتصال بالمنصات
const ALLOW = ["https://abd6996.github.io"];
export default {
  async fetch(req) {
    const o = req.headers.get("Origin") || "";
    const cors = {"Access-Control-Allow-Origin": o, "Access-Control-Allow-Headers": "x-h,content-type", "Access-Control-Expose-Headers": "x-final-url,content-length,content-range,content-type", "Vary": "Origin"};
    if (req.method === "OPTIONS") return new Response(null, {headers: cors});
    if (!ALLOW.includes(o)) return new Response("forbidden", {status: 403});
    const u = new URL(req.url).searchParams.get("u") || "";
    if (!/^https?:\/\//i.test(u)) return new Response("bad url", {status: 400, headers: cors});
    let h = {};
    try { h = JSON.parse(req.headers.get("x-h") || "{}"); } catch (e) {}
    try {
      const r = await fetch(u, {method: req.method, headers: h, body: ["GET", "HEAD"].includes(req.method) ? undefined : await req.text(), redirect: "follow"});
      const rh = new Headers(r.headers);
      rh.delete("set-cookie");
      for (const k in cors) rh.set(k, cors[k]);
      rh.set("x-final-url", r.url);
      return new Response(r.body, {status: r.status, headers: rh});
    } catch (e) { return new Response("fetch failed", {status: 502, headers: cors}); }
  }
};
