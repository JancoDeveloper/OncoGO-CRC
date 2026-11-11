var Comparador = (function () {
  function escapar(texto) {
    return String(texto == null ? "" : texto)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  }
  function conjunto(lista) {
    var sale = {};
    (lista || []).forEach(function (item) { sale[item.id] = item; });
    return sale;
  }
  function jaccard(a, b) {
    var claves = {};
    Object.keys(a).forEach(function (k) { claves[k] = true; });
    Object.keys(b).forEach(function (k) { claves[k] = true; });
    var union = Object.keys(claves);
    if (!union.length) return 0;
    var ambos = union.filter(function (k) { return a[k] && b[k]; }).length;
    return ambos / union.length;
  }
  function fila(celdas) {
    return "<tr>" + celdas.map(function (celda) { return "<td>" + celda + "</td>"; }).join("") + "</tr>";
  }
  function bloque(id, titulo, cuerpo) {
    return '<h2 id="' + id + '">' + escapar(titulo) + "</h2>" + cuerpo;
  }
  function pintar(ids, catalogo) {
    var proteínas = ids.map(function (id) { return catalogo[id]; }).filter(Boolean);
    if (proteínas.length < 2) return "<p>Hacen falta de 2 a 4 accesos que estén en la celda.</p>";
    var resumen = "<table><thead><tr><th>Campo</th>" + proteínas.map(function (p) {
      return "<th>" + escapar(p.gen || p.id) + "</th>";
    }).join("") + "</tr></thead><tbody>";
    [["UniProt", "revisada"], ["Longitud", "longitud"], ["Masa", "masa"], ["Clase CP frente a CN", "clase"], ["log2FC CP frente a CN", "log2fc_cn"], ["Paso", "paso"]].forEach(function (par) {
      resumen += fila(["<th>" + escapar(par[0]) + "</th>"].concat(proteínas.map(function (p) {
        return escapar(p[par[1]] == null ? "—" : p[par[1]]);
      })));
    });
    resumen += "</tbody></table>";
    var go = "";
    for (var i = 0; i < proteínas.length; i++) {
      for (var j = i + 1; j < proteínas.length; j++) {
        var a = conjunto(proteínas[i].go);
        var b = conjunto(proteínas[j].go);
        var ap = conjunto(proteínas[i].go_propagado);
        var bp = conjunto(proteínas[j].go_propagado);
        go += "<p>" + escapar(proteínas[i].id) + " y " + escapar(proteínas[j].id)
          + ". Jaccard directo " + jaccard(a, b).toFixed(3)
          + ". Jaccard propagado " + jaccard(ap, bp).toFixed(3) + ".</p>";
      }
    }
    var secuencia = "";
    for (i = 0; i < proteínas.length; i++) {
      for (j = i + 1; j < proteínas.length; j++) {
        var aln = Alineamiento.alinear(proteínas[i].secuencia, proteínas[j].secuencia);
        var identidad = aln.columnas ? (100 * aln.identidades / aln.columnas) : 0;
        secuencia += "<h3>" + escapar(proteínas[i].id) + " / " + escapar(proteínas[j].id) + "</h3>"
          + "<p>Puntuación " + aln.puntuacion + ". Identidad " + identidad.toFixed(1) + " %. "
          + "Similitud " + (aln.columnas ? (100 * aln.similares / aln.columnas) : 0).toFixed(1) + " %.</p>"
          + '<p class="secuencia">' + escapar(aln.arriba) + "\n" + escapar(aln.abajo) + "</p>";
      }
    }
    var dominios = "<table><thead><tr><th>Acceso</th><th>InterPro</th><th>Pfam</th><th>Reactome</th></tr></thead><tbody>";
    proteínas.forEach(function (p) {
      dominios += fila([
        escapar(p.id),
        escapar((p.interpro || []).join("; ")),
        escapar((p.pfam || []).join("; ")),
        escapar((p.reactome || []).join("; "))
      ]);
    });
    dominios += "</tbody></table>";
    var estructura = "<table><thead><tr><th>Acceso</th><th>pLDDT medio</th></tr></thead><tbody>";
    var modelos = [];
    proteínas.forEach(function (p) {
      estructura += fila([escapar(p.id), escapar(p.plddt == null ? "sin modelo" : p.plddt)]);
      if (p.modelo) modelos.push(p.modelo);
    });
    estructura += "</tbody></table>";
    if (modelos.length >= 2) {
      estructura += '<p><button type="button" data-ver-estructura="visor-par" data-modelos="'
        + escapar(modelos.slice(0, 2).join(","))
        + '" data-script="' + escapar(catalogo._molstar_js)
        + '" data-estilo="' + escapar(catalogo._molstar_css)
        + '">Ver estructuras</button></p><div id="visor-par"></div>';
    } else {
      estructura += "<p>No hay dos modelos alojados para abrirlos juntos. No se superponen estructuras.</p>";
    }
    return bloque("vista-resumen", "Resumen", resumen)
      + bloque("vista-go", "GO", go)
      + bloque("vista-secuencia", "Secuencia", secuencia)
      + bloque("vista-dominios", "Dominios", dominios)
      + bloque("vista-estructura", "Estructura", estructura);
  }
  function leerConsulta() {
    var params = new URLSearchParams(window.location.search);
    return (params.get("p") || "").split(",").map(function (id) { return id.trim(); }).filter(Boolean);
  }
  function iniciar() {
    var caja = document.getElementById("comparacion");
    if (!caja) return;
    var pedido = leerConsulta();
    var actual = (caja.getAttribute("data-par") || "").split(",");
    if (pedido.length === 1) {
      var campo = document.getElementById("p");
      if (campo) campo.value = pedido[0] + ",";
      return;
    }
    if (pedido.length < 2 || pedido.join(",") === actual.join(",")) return;
    fetch("/fuentes/comparador.json").then(function (respuesta) { return respuesta.json(); }).then(function (datos) {
      datos.proteinas._molstar_js = datos.molstar_js;
      datos.proteinas._molstar_css = datos.molstar_css;
      caja.innerHTML = pintar(pedido.slice(0, 4), datos.proteinas);
    });
  }
  return { iniciar: iniciar };
})();
Comparador.iniciar();
