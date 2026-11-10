
document.querySelectorAll("details").forEach((nodo) => {
  nodo.addEventListener("toggle", () => nodo.setAttribute("data-abierto", nodo.open ? "si" : "no"));
});
