(function () {
  var raiz = document.documentElement;
  var guardado = null;
  try { guardado = localStorage.getItem("tema"); } catch (error) { guardado = null; }
  if (guardado === "oscuro" || guardado === "claro") raiz.setAttribute("data-tema", guardado);
  var boton = document.getElementById("tema");
  if (boton) {
    boton.addEventListener("click", function () {
      var siguiente = raiz.getAttribute("data-tema") === "oscuro" ? "claro" : "oscuro";
      raiz.setAttribute("data-tema", siguiente);
      try { localStorage.setItem("tema", siguiente); } catch (error) { return; }
    });
  }
  document.querySelectorAll("img[data-respaldo]").forEach(function (img) {
    img.addEventListener("error", function () {
      if (img.getAttribute("data-cambio") === "1") return;
      img.setAttribute("data-cambio", "1");
      img.src = img.getAttribute("data-respaldo");
    });
  });
  var path = location.pathname.replace(/\/$/, "") || "/";
  document.querySelectorAll(".menu a").forEach(function (enlace) {
    var href = enlace.getAttribute("href");
    var activo = href === "/" ? path === "/" : path === href || path.indexOf(href + "/") === 0;
    if (activo) enlace.setAttribute("aria-current", "page");
  });
  var caja = document.getElementById("q");
  if (caja && document.querySelector("[data-buscar]")) {
    var aplicar = function () {
      var consulta = caja.value.toLowerCase();
      document.querySelectorAll("[data-buscar]").forEach(function (nodo) {
        if (nodo.closest && nodo.closest(".explorador")) return;
        var texto = nodo.getAttribute("data-buscar").toLowerCase();
        nodo.hidden = consulta.length > 0 && texto.indexOf(consulta) === -1;
      });
    };
    caja.addEventListener("input", aplicar);
    if (caja.value) aplicar();
  }
  var secciones = Array.prototype.slice.call(document.querySelectorAll("main h2[id]"));
  var enlaces = Array.prototype.slice.call(document.querySelectorAll(".indice a"));
  if (secciones.length && enlaces.length && "IntersectionObserver" in window) {
    var observador = new IntersectionObserver(function (entradas) {
      entradas.forEach(function (entrada) {
        if (!entrada.isIntersecting) return;
        enlaces.forEach(function (enlace) {
          var activo = enlace.getAttribute("href") === "#" + entrada.target.id;
          if (activo) enlace.setAttribute("aria-current", "true");
          else enlace.removeAttribute("aria-current");
        });
      });
    }, { rootMargin: "-20% 0px -70% 0px" });
    secciones.forEach(function (seccion) { observador.observe(seccion); });
  }
  document.querySelectorAll("[data-filtro]").forEach(function (boton) {
    boton.addEventListener("click", function () {
      var filtro = boton.getAttribute("data-filtro");
      var valor = boton.getAttribute("data-valor");
      document.querySelectorAll(".celda-multi").forEach(function (celda) {
        var ocultar = valor !== "todas" && celda.getAttribute(filtro) !== valor;
        celda.style.display = ocultar ? "none" : "";
      });
    });
  });
  document.body.addEventListener("click", function (evento) {
    var boton = evento.target.closest ? evento.target.closest("[data-ver-estructura]") : null;
    if (!boton) return;
    var destino = document.getElementById(boton.getAttribute("data-ver-estructura"));
    if (!destino || destino.getAttribute("data-listo") === "1") return;
    destino.setAttribute("data-listo", "1");
    var estilo = document.createElement("link");
    estilo.rel = "stylesheet";
    estilo.href = boton.getAttribute("data-estilo");
    document.head.appendChild(estilo);
    var cargador = document.createElement("script");
    cargador.src = boton.getAttribute("data-script");
    cargador.onload = function () {
      if (!window.molstar || !molstar.Viewer) {
        destino.textContent = "El visor no está disponible en esta página. La cifra de pLDDT sigue en el texto.";
        return;
      }
      var archivos = (boton.getAttribute("data-modelos") || "").split(",").filter(Boolean);
      archivos.forEach(function (archivo) {
        var caja = document.createElement("div");
        caja.className = "visor";
        destino.appendChild(caja);
        molstar.Viewer.create(caja, {
          layoutIsExpanded: false,
          layoutShowControls: false,
          layoutShowSequence: true,
          layoutShowLog: false,
          layoutShowLeftPanel: false,
          viewportShowExpand: false,
          viewportShowSettings: false,
          viewportShowSelectionMode: false,
          viewportShowAnimation: false,
          viewportShowTrajectoryControls: false
        }).then(function (visor) {
          return visor.loadStructureFromUrl(archivo, "mmcif", true).then(function () {
            var plugin = visor.plugin;
            if (plugin && plugin.canvas3d && plugin.canvas3d.handle && plugin.canvas3d.handle.resize) {
              plugin.canvas3d.handle.resize();
            }
          });
        }).catch(function () {
          caja.textContent = "El modelo no se abrió en este navegador. La cifra de pLDDT sigue en el texto.";
        });
      });
    };
    cargador.onerror = function () {
      destino.textContent = "No se pudo cargar el visor. La cifra de pLDDT sigue en el texto.";
    };
    document.body.appendChild(cargador);
  });
  document.querySelectorAll("[data-segmento]").forEach(function (boton) {
    boton.addEventListener("click", function () {
      var id = boton.getAttribute("aria-controls");
      document.querySelectorAll("[data-segmento]").forEach(function (otro) {
        otro.setAttribute("aria-pressed", otro === boton ? "true" : "false");
      });
      document.querySelectorAll(".segmento").forEach(function (panel) {
        panel.hidden = panel.id !== id;
      });
    });
  });
  document.querySelectorAll(".embudo, .matrices").forEach(function (bloque) {
    if (window.matchMedia("(prefers-reduced-motion: reduce)").matches) return;
    bloque.classList.add("animar");
    window.requestAnimationFrame(function () { bloque.classList.add("listo"); });
  });
  var lista = document.querySelector(".explorador");
  if (lista) {
    var items = Array.prototype.slice.call(lista.querySelectorAll("li"));
    var filtro = document.getElementById("filtro-proteinas");
    var solo = document.getElementById("solo-paso-5");
    var pagina = 0;
    var tam = 24;
    var estado = document.createElement("p");
    var controles = document.createElement("p");
    controles.className = "paginas";
    var anterior = document.createElement("button");
    var siguiente = document.createElement("button");
    anterior.type = "button";
    siguiente.type = "button";
    anterior.textContent = "Anterior";
    siguiente.textContent = "Siguiente";
    controles.appendChild(anterior);
    controles.appendChild(siguiente);
    lista.after(estado, controles);
    function visibles() {
      var q = filtro && filtro.value ? filtro.value.trim().toLowerCase() : "";
      var paso = solo && solo.checked;
      return items.filter(function (li) {
        var texto = (li.getAttribute("data-buscar") || "").toLowerCase();
        if (paso && texto.indexOf("candidata") === -1) return false;
        return !q || texto.indexOf(q) !== -1;
      });
    }
    function pintar() {
      var listaVisible = visibles();
      var paginas = Math.max(1, Math.ceil(listaVisible.length / tam));
      if (pagina >= paginas) pagina = paginas - 1;
      if (pagina < 0) pagina = 0;
      items.forEach(function (li) { li.hidden = true; });
      listaVisible.slice(pagina * tam, pagina * tam + tam).forEach(function (li) { li.hidden = false; });
      var desde = listaVisible.length ? pagina * tam + 1 : 0;
      var hasta = Math.min(listaVisible.length, (pagina + 1) * tam);
      estado.textContent = desde + "–" + hasta + " de " + listaVisible.length;
      anterior.disabled = pagina === 0;
      siguiente.disabled = pagina >= paginas - 1;
    }
    anterior.addEventListener("click", function () { pagina -= 1; pintar(); });
    siguiente.addEventListener("click", function () { pagina += 1; pintar(); });
    if (filtro) filtro.addEventListener("input", function () { pagina = 0; pintar(); });
    if (solo) solo.addEventListener("change", function () { pagina = 0; pintar(); });
    pintar();
  }
  document.querySelectorAll(".lienzo-svg [data-id]").forEach(function (nodo) {
    nodo.addEventListener("click", function () {
      var svg = nodo.closest("svg");
      svg.querySelectorAll(".resaltado").forEach(function (el) { el.classList.remove("resaltado"); });
      nodo.classList.add("resaltado");
      svg.querySelectorAll("line").forEach(function (linea) {
        var id = nodo.getAttribute("data-id");
        if (linea.getAttribute("data-origen") === id || linea.getAttribute("data-destino") === id) {
          linea.classList.add("resaltado");
        }
      });
    });
  });
})();
