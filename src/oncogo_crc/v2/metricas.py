"""Fracciones con numerador y denominador obligatorios (corrección C3).

En v1 el certificado publicaba 0,3285 y 0,1695 como «la misma» fracción de
detectabilidad con denominadores distintos: 688 proteínas en un caso y 702
**filas** del cruce en el otro, porque una proteína de CPTAC puede mapear a
varios UniProt. Un número suelto no deja ver eso.

Aquí una fracción no se puede construir sin decir qué cuenta arriba, qué cuenta
abajo y en qué unidad. El valor se deriva; nunca se pasa a mano.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from scipy.stats import binomtest

UNIDADES_VALIDAS = frozenset({"proteinas_unicas", "terminos", "permutaciones", "pacientes", "pares"})


@dataclass(frozen=True)
class Fraccion:
    """Una proporción que sabe de dónde salió.

    ``cobertura`` es opcional y sirve para R1: dice qué parte del universo de
    interés entra en el denominador. Sin ella, una concordancia medida solo en
    las proteínas que la réplica alcanza a ver se lee como si fuera de todas.
    """

    nombre: str
    numerador: int
    denominador: int
    unidad_denominador: str
    universo: int | None = None
    universo_nombre: str = ""
    nota: str = ""
    extra: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.unidad_denominador not in UNIDADES_VALIDAS:
            raise ValueError(
                f"unidad {self.unidad_denominador!r} no válida; use una de {sorted(UNIDADES_VALIDAS)}"
            )
        if self.denominador < 0 or self.numerador < 0:
            raise ValueError(f"{self.nombre}: numerador y denominador no pueden ser negativos")
        if self.numerador > self.denominador:
            raise ValueError(
                f"{self.nombre}: numerador {self.numerador} > denominador {self.denominador}"
            )

    @property
    def valor(self) -> float | None:
        if self.denominador == 0:
            return None
        return self.numerador / self.denominador

    @property
    def cobertura(self) -> float | None:
        """Denominador sobre universo. Para R1: 119/688 = 0,173."""
        if self.universo in (None, 0):
            return None
        return self.denominador / float(self.universo)

    def intervalo(self, confianza: float = 0.95) -> tuple[float | None, float | None]:
        """Clopper–Pearson exacto. Con denominador cero no hay intervalo."""
        if self.denominador == 0:
            return (None, None)
        test = binomtest(self.numerador, self.denominador, p=0.5)
        ci = test.proportion_ci(confidence_level=confianza, method="exact")
        return (float(ci.low), float(ci.high))

    def p_vs_azar(self) -> float | None:
        if self.denominador == 0:
            return None
        return float(binomtest(self.numerador, self.denominador, p=0.5).pvalue)

    def to_dict(self) -> dict[str, Any]:
        low, high = self.intervalo()
        payload: dict[str, Any] = {
            "nombre": self.nombre,
            "valor": self.valor,
            "numerador": self.numerador,
            "denominador": self.denominador,
            "unidad_denominador": self.unidad_denominador,
            "ci_low": low,
            "ci_high": high,
            "p_vs_azar": self.p_vs_azar(),
        }
        if self.universo is not None:
            payload["universo"] = self.universo
            payload["universo_nombre"] = self.universo_nombre
            payload["cobertura"] = self.cobertura
        if self.nota:
            payload["nota"] = self.nota
        payload.update(self.extra)
        return payload


def p_empirico(b: int, n_permutaciones: int) -> dict[str, Any]:
    """p de permutación con el piso declarado (corrección C1).

    ``b`` es cuántas permutaciones igualaron o superaron lo observado. El p es
    (b+1)/(B+1), que nunca baja de 1/(B+1). Cuando b = 0 el valor devuelto es ese
    piso y ``es_tope`` lo dice, para que el reporte no lo presente como una
    estimación exacta.
    """
    if n_permutaciones <= 0:
        raise ValueError("n_permutaciones debe ser positivo")
    frac = Fraccion(
        nombre="p_empirico",
        numerador=b + 1,
        denominador=n_permutaciones + 1,
        unidad_denominador="permutaciones",
    )
    low, high = frac.intervalo()
    return {
        "p": frac.valor,
        "b": int(b),
        "n_permutaciones": int(n_permutaciones),
        "p_minimo_alcanzable": 1.0 / (n_permutaciones + 1),
        "es_tope": b == 0,
        "ci_low": low,
        "ci_high": high,
        "formula": "(b + 1) / (B + 1)",
    }
