"""El mismo formato numérico para la tesis, la web y las figuras."""

from oncogo_crc.formato import formatear


def test_auc_con_tres_decimales_y_la_celda_completa():
    assert formatear(0.88025, "auc") == "0,880"
    assert formatear(0.9151250000000001, "auc") == "0,915"
    assert formatear(0.64675, "auc") == "0,647"
    assert formatear(0.88025, "auc_completo") == "0,88025"


def test_p_con_dos_cifras_y_cientifica_por_debajo_de_un_milesimo():
    assert formatear(0.11888111888111888, "p") == "0,12"
    assert formatear(0.00099, "p") == "0,00099"
    assert formatear(2.1733818203621972e-08, "p") == "2,2 × 10⁻⁸"


def test_log2fc_rho_enteros_y_estado():
    assert formatear(3.027, "log2fc") == "3,03"
    assert formatear(-0.77915, "rho") == "−0,78"
    assert formatear(2934, "entero") == "2\u202f934"
    assert formatear("no_evaluable", "estado") == "No evaluable"
