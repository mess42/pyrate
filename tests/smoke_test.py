#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2018
               by     Moritz Esslinger moritz.esslinger@web.de
               and    Johannes Hartung j.hartung@gmx.net
               and    Uwe Lippmann  uwe.lippmann@web.de

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""

from matplotlib.testing.decorators import image_comparison
from matplotlib import pyplot

# we use the matplotlib's do nothing backend for testing
# matplotlib/lib/matplotlib/backends/backend_template.py
pyplot.switch_backend('Template')

def test_smoke_anisotropic_doublet():
    """Smoke test based on demo_anisotropic_doublet.py."""
    import demos.demo_anisotropic_doublet
    assert True
def test_smoke_anisotropic_mirror():
    """Smoke test based on demo_anisotropic_mirror.py."""
    import demos.demo_anisotropic_mirror
    assert True
def test_smoke_anisotropic_mirrors():
    """Smoke test based on demo_anisotropic_mirrors.py."""
    import demos.demo_anisotropic_mirrors
    assert True
def test_smoke_anisotropic_ord_eo():
    """Smoke test based on demo_anisotropic_ord_eo.py."""
    import demos.demo_anisotropic_ord_eo
    assert True
def test_smoke_asphere():
    """Smoke test based on demo_asphere.py."""
    import demos.demo_asphere
    assert True
# TODO: currently not working, due to missing refractingindex database
# def test_smoke_benchmark():
#    """Smoke test based on demo_benchmark.py."""
#    import demos.demo_benchmark
#    assert True
# TODO: currently not working, due to missing refractingindex database
# def test_smoke_doublegauss():
#    """Smoke test based on demo_doublegauss.py."""
#    import demos.demo_doublegauss
#    assert True
@image_comparison(baseline_images=['doublet'], extensions=['png'])
def test_smoke_doublet():
    """Smoke test based on demo_doublet.py."""
    import demos.demo_doublet
    assert True
def test_smoke_grin():
    """Smoke test based on demo_grin.py."""
    import demos.demo_grin
    assert True
def test_smoke_hud():
    """Smoke test based on demo_hud.py."""
    import demos.demo_hud
    assert True
# TODO: currently not working, outdated demo script ?
# def test_smoke_loadsave():
#    """Smoke test based on demo_loadsave.py."""
#    import demos.demo_loadsave
#    assert True
def test_smoke_mirrors():
    """Smoke test based on demo_mirrors.py."""
    import demos.demo_mirrors
    assert True
def test_smoke_optimize():
    """Smoke test based on demo_optimize.py."""
    import demos.demo_optimize
    assert True
def test_smoke_prism():
    """Smoke test based on demo_prism.py."""
    import demos.demo_prism
    assert True
def test_smoke_rainbow():
    """Smoke test based on demo_rainbow.py."""
    import demos.demo_rainbow
    assert True
def test_smoke_zmx():
    """Smoke test based on demo_zmx.py."""
    import sys
    sys.argv = [sys.argv[0], 'tests/lenssystem.ZMX', '10', '2']
    import demos.demo_zmx
    assert True
