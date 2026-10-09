import numpy as np
import dexter as dex


def test_lar_stationary_curve(lar_machine: dex.Machine):
    curve = dex.StationaryCurve(lar_machine)
    assert curve.flux_kind == "Toroidal"
    segments = curve.segments()
    assert len(segments) == 2
    for seg in segments:
        flux = seg.flux
        theta = seg.theta
        assert isinstance(flux, np.ndarray)
        assert isinstance(theta, np.ndarray)
        assert len(flux) == len(theta) == len(seg)
