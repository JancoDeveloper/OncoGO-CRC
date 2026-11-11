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
  var caja = document.getElementById("q");
  if (caja && document.querySelector("[data-buscar]")) {
    var aplicar = function () {
      var consulta = caja.value.toLowerCase();
      document.querySelectorAll("[data-buscar]").forEach(function (nodo) {
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
          return visor.loadStructureFromUrl(archivo, "bcif", true);
        });
      });
    };
    document.body.appendChild(cargador);
  });
  document.querySelectorAll(".lienzo-svg a[data-id]").forEach(function (nodo) {
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
