import numpy as np
import pandas as pd
import pytest

from deprem.mekanik import termal


def test_jeoterm_ve_ters_cozum():
    q = termal.isi_akisi_derinlikten(14.0)
    assert 0.05 < q < 0.12
    assert abs(termal.sicaklik(14.0, q) - termal.T_KIRILGAN) < 1e-6
    # daha sıcak (yüksek ısı akısı) -> daha ince kırılgan katman
    assert termal.isi_akisi_derinlikten(10.0) > termal.isi_akisi_derinlikten(16.0)


def test_surtunme_rejimi():
    ab = termal.a_eksi_b([20, 200, 500])
    assert ab[0] > 0 and ab[1] < 0 and ab[2] > 0


def _kucuk_fay():
    cutde = pytest.importorskip("cutde")  # noqa: F841
    from deprem.mekanik import geometri as g
    iz = [{"ad": "test", "coords": np.array([[29.0, 40.7], [29.4, 40.7]])}]
    c = g.hucreler(iz, z_kirilgan=12.0, z_sunek=20.0)
    Kt, Kn = g.etkilesim_matrisi(c)
    return c, Kt, Kn


def test_etkilesim_matrisi_fiziksel():
    c, Kt, Kn = _kucuk_fay()
    i = c.index[(c.z_km == 5.5) & (c.x_km.between(16, 17))][0]
    # kendi kayması gerilmeyi düşürür, komşuyu yükler
    assert Kt[i, i] < -2e7
    j = c.index[(c.z_km == 5.5) & (c.x_km.between(17, 18))][0]
    assert Kt[i, j] > 0
    # düzlemsel düşey fayda kayma normal gerilme üretmez (yaklaşık)
    assert np.abs(Kn[i]).max() < 0.05 * abs(Kt[i, i])
    # tüm şerit düzgün kayarsa gerilme, kendi kendine etkiden çok küçüktür
    assert abs(Kt[i].sum()) < 0.1 * abs(Kt[i, i])


def test_simulator_deprem_uretir_ve_kayma_birikir():
    from deprem.mekanik.simulator import Model, Parametreler
    c, Kt, Kn = _kucuk_fay()
    m = Model(c, Kt, Kn, Parametreler(dtau0=3e6, seed=1))
    m.ilerle(400)
    cat = m.katalog()
    assert len(cat) > 5
    assert cat.mw.max() > 5.5
    # uzun dönemde kilitli hücrelerin kayması levha hızına yaklaşır (deprem + krip)
    oran = m.slip[m.vw].mean() / (m.v_pl[m.vw].mean() * 400)
    assert 0.6 < oran < 1.3
