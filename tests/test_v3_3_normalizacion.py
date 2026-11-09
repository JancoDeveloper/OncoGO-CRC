import numpy as np
import pandas as pd

from oncogo_crc.v3.etiquetas import centrar_por_mediana


def test_centrar_recupera_un_log2fc_cerca_de_cero():
    mat = pd.DataFrame(np.full((120, 4), 5.0), columns=["P1", "P2", "N1", "N2"])
    mat[["P1", "P2"]] = mat[["P1", "P2"]] + 1.0
    centrada, _desplazamiento, regla = centrar_por_mediana(mat)
    assert regla == "100 %"
    diferencia = (centrada["P1"] - centrada["N1"]).median()
    assert abs(float(diferencia)) < 1e-9
