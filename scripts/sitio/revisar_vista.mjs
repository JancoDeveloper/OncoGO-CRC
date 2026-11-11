#!/usr/bin/env node
/* Revisa el sitio ya servido. No inventa páginas: lee los HTML de sitio/. */

const fs = require("fs");
const path = require("path");
const { chromium } = require("playwright");
const { AxeBuilder } = require("@axe-core/playwright");

const raiz = process.cwd();
const base = process.env.BASE || "http://127.0.0.1:5002";
const capturas = path.join(raiz, "docs", "sitio_v4_capturas");

function rutas() {
  const sitio = path.join(raiz, "sitio");
  const salida = [];
  function caminar(dir) {
    for (const nombre of fs.readdirSync(dir)) {
      if (nombre === "archivo") continue;
      const completo = path.join(dir, nombre);
      const stat = fs.statSync(completo);
      if (stat.isDirectory()) caminar(completo);
      else if (nombre.endsWith(".html")) {
        const rel = path.relative(sitio, completo).split(path.sep).join("/");
        const partes = rel === "index.html"
          ? []
          : rel.replace(/\/index\.html$/, "").replace(/\.html$/, "").split("/");
        const url = "/" + partes.map((parte) => encodeURIComponent(parte).replace(/%3B/g, ";")).join("/");
        salida.push(url);
      }
    }
  }
  caminar(sitio);
  return salida.sort();
}

async function axe(pagina, url) {
  await pagina.goto(base + url, { waitUntil: "domcontentloaded" });
  const resultado = await new AxeBuilder({ page: pagina }).analyze();
  return resultado.violations.filter((item) => item.impact === "serious" || item.impact === "critical");
}

async function foto(pagina, nombre, url, ancho, tema, preparar) {
  await pagina.setViewportSize({ width: ancho, height: ancho === 390 ? 844 : 900 });
  await pagina.goto(base + url, { waitUntil: "domcontentloaded" });
  await pagina.evaluate((valor) => {
    document.documentElement.setAttribute("data-tema", valor);
  }, tema);
  if (preparar) await preparar(pagina);
  const destino = path.join(capturas, `${nombre}-${ancho}-${tema}.png`);
  await pagina.screenshot({ path: destino });
}

async function abrirEstructura(pagina) {
  const boton = pagina.locator("[data-ver-estructura]").first();
  await boton.click();
  await pagina.locator("#visor-P02749 canvas").first().waitFor({ timeout: 30000 });
  await pagina.waitForTimeout(2500);
  await pagina.waitForFunction(() => {
    const canvas = document.querySelector("#visor-P02749 canvas");
    const gl = canvas && (canvas.getContext("webgl2") || canvas.getContext("webgl"));
    if (!gl) return false;
    let oscuros = 0;
    for (let x = 0; x < canvas.width; x += 16) {
      for (let y = 0; y < canvas.height; y += 16) {
        const buf = new Uint8Array(4);
        gl.readPixels(x, y, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, buf);
        if (buf[0] < 230 || buf[1] < 230 || buf[2] < 230) oscuros += 1;
      }
    }
    return oscuros > 3;
  }, null, { timeout: 20000 });
  await pagina.evaluate(() => {
    const visor = document.querySelector("#visor-P02749 .visor");
    const canvas = visor.querySelector("canvas");
    const gl = canvas.getContext("webgl2") || canvas.getContext("webgl");
    const w = canvas.width;
    const h = canvas.height;
    const pixels = new Uint8Array(w * h * 4);
    gl.readPixels(0, 0, w, h, gl.RGBA, gl.UNSIGNED_BYTE, pixels);
    const copia = document.createElement("canvas");
    copia.width = w;
    copia.height = h;
    copia.style.width = "100%";
    copia.style.height = "100%";
    const ctx = copia.getContext("2d");
    const img = ctx.createImageData(w, h);
    for (let y = 0; y < h; y += 1) {
      img.data.set(pixels.subarray((h - 1 - y) * w * 4, (h - y) * w * 4), y * w * 4);
    }
    ctx.putImageData(img, 0, 0);
    visor.replaceChildren(copia);
    visor.scrollIntoView({ block: "center" });
  });
}

async function main() {
  fs.mkdirSync(capturas, { recursive: true });
  const navegador = await chromium.launch({
    headless: true,
    args: ["--use-gl=angle", "--use-angle=swiftshader", "--enable-webgl", "--ignore-gpu-blocklist"],
  });
  const contexto = await navegador.newContext();
  const pagina = await contexto.newPage();
  if (process.env.SOLO_ESTRUCTURA === "1") {
    for (const ancho of [1440, 390]) {
      for (const tema of ["claro", "oscuro"]) {
        await foto(pagina, "p02749-estructura", "/proteinas/P02749", ancho, tema, abrirEstructura);
        console.log("foto", "p02749-estructura", ancho, tema);
      }
    }
    await navegador.close();
    return;
  }
  const fallos = [];
  const lista = rutas();
  for (let i = 0; i < lista.length; i += 1) {
    const url = lista[i];
    const graves = await axe(pagina, url);
    if (graves.length) {
      fallos.push({
        url,
        reglas: graves.map((item) => ({
          id: item.id,
          impact: item.impact,
          nodos: item.nodes.length,
          ayuda: item.help,
        })),
      });
      console.log("axe", url, graves.map((item) => item.id).join(","));
    } else if (i % 40 === 0) {
      console.log("va", i + 1, "/", lista.length);
    }
  }
  const informe = path.join(raiz, "docs", "sitio_v4_axe.json");
  fs.writeFileSync(informe, JSON.stringify({ base, paginas: lista.length, fallos }, null, 2) + "\n");
  if (fallos.length) {
    await navegador.close();
    console.log("fallos", fallos.length);
    process.exit(1);
  }
  const vistas = [
    ["portada", "/"],
    ["contexto", "/investigacion/contexto"],
    ["historia", "/investigacion/historia"],
    ["datos", "/datos"],
    ["bases-de-datos", "/datos/bases-de-datos"],
    ["canalizacion", "/metodos/canalizacion"],
    ["expresion", "/metodos/expresion-diferencial"],
    ["suero", "/resultados/suero"],
    ["discusion", "/discusion"],
    ["glosario", "/glosario"],
    ["bibliografia", "/bibliografia"],
    ["proteinas", "/proteinas"],
    ["comparador", "/proteinas/comparar?p=P02749,P05155"],
    ["candidatas", "/proteinas/candidatas"],
  ];
  for (const [nombre, url] of vistas) {
    for (const ancho of [1440, 390]) {
      for (const tema of ["claro", "oscuro"]) {
        await foto(pagina, nombre, url, ancho, tema);
        console.log("foto", nombre, ancho, tema);
      }
    }
  }
  for (const ancho of [1440, 390]) {
    for (const tema of ["claro", "oscuro"]) {
      await foto(pagina, "p02749-estructura", "/proteinas/P02749", ancho, tema, abrirEstructura);
      console.log("foto", "p02749-estructura", ancho, tema);
    }
  }
  await navegador.close();
  console.log("listo", lista.length);
}

main().catch((error) => {
  console.error(error);
  process.exit(1);
});
