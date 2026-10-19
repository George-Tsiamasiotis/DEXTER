import dexter as dex


def test_stationary_curve_creation_analytical(lar_machine: dex.Machine):
    curve = dex.StationaryCurve(lar_machine)
    segments = curve.segments
    assert len(segments) == len(curve) == 2


def test_stationary_curve_creation_numerical(nc_machine: dex.Machine):
    curve = dex.StationaryCurve(nc_machine)
    segments = curve.segments
    assert len(segments) == len(curve) == 2
