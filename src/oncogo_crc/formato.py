"""Un solo formato numérico para la tesis, la web y las figuras.

Los JSON guardan el valor completo. Esta función solo decide cómo se ve.
"""

from __future__ import annotations

import math

_ESTADOS = {
    "sostiene": "Sostiene",
    "no_sostiene": "No sostiene",
    "refuta": "Refuta",
    "no_evaluable": "No evaluable",
    "exploratorio": "Exploratorio",
    "descriptiva": "Descriptiva",
    "calculado": "Calculado",
    "pasa": "Pasa",
    "falla": "Falla",
    "publicable": "Publicable",
    "no_disponible": "No disponible",
    "no_detectado": "No detectado",
    "no_medida": "No medida",
}
_SUPER = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
_MILES = "\u202f"


def formatear(valor: object, tipo: str) -> str:
    if valor is None or (isinstance(valor, float) and math.isnan(valor)):
        return "no medido"
    if tipo == "estado":
        texto = str(valor)
        return _ESTADOS.get(texto, texto.replace("_", " ").capitalize())
    numero = float(valor)
    if tipo == "auc":
        return _coma(f"{numero:.3f}")
    if tipo == "auc_completo":
        return _coma(f"{numero:.5f}".rstrip("0").rstrip("."))
    if tipo == "p":
        return _p(numero)
    if tipo in {"log2fc", "beta", "rho"}:
        return _coma(f"{abs(numero):.2f}", negativo=numero < 0)
    if tipo == "porcentaje":
        decimales = 1 if abs(numero) < 10 and not math.isclose(numero, round(numero)) else 0
        return _coma(f"{numero:.{decimales}f}") + " %"
    if tipo == "entero":
        return _miles(int(round(numero)))
    raise ValueError(tipo)


def _coma(texto: str, negativo: bool = False) -> str:
    if texto.startswith("-"):
        negativo = True
        texto = texto[1:]
    texto = texto.replace(".", ",")
    return ("−" if negativo else "") + texto


def _miles(valor: int) -> str:
    signo = "−" if valor < 0 else ""
    digitos = f"{abs(valor):,}".replace(",", _MILES)
    return signo + digitos


def _p(valor: float) -> str:
    if valor == 0:
        return "0"
    if valor < 0:
        return "−" + _p(-valor)
    if valor < 1e-4:
        mantisa, exponente = f"{valor:.1e}".split("e")
        return f"{mantisa.replace('.', ',')} × 10{str(int(exponente)).translate(_SUPER)}"
    texto = f"{valor:.2g}"
    if "e" in texto:
        mantisa, exponente = texto.split("e")
        return f"{mantisa.replace('.', ',')} × 10{exponente.lstrip('+').translate(_SUPER)}"
    return _coma(texto)
