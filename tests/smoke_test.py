"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2017 Moritz Esslinger <moritz.esslinger@web.de>
               and Johannes Hartung <j.hartung@gmx.net>
               and     Uwe Lippmann <uwe.lippmann@web.de>
               and    Thomas Heinze <t.heinze@fn.de>

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
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA.
"""

from matplotlib.testing.decorators import image_comparison
from matplotlib import pyplot

__author__ = 'Thomas Heinze'

# we use the matplotlib's do nothing backend for testing
# matplotlib/lib/matplotlib/backends/backend_template.py
pyplot.switch_backend('Template')

@image_comparison(baseline_images=['doublet'], extensions=['png'])
def test_smoke_doublet():
    """Smoke test based on demo_doublet.py."""
    import demo_doublet
    assert True

def test_smoke_benchmark():
    """Smoke test based on demo_benchmark.py."""
    import demo_benchmark
    assert True
def test_smoke_mirror():
    """Smoke test based on demo_mirrors.py."""
    import demo_mirrors
    assert True
def test_smoke_optimize():
    """Smoke test based on demo_optimize.py."""
    import demo_optimize
    assert True
def test_smoke_grin():
    """Smoke test based on demo_grin.py."""
    import demo_grin
    assert True
def test_smoke_prism():
    """Smoke test based on demo_prism.py."""
    import demo_prism
    assert True
def test_smoke_rainbow():
    """Smoke test based on demo_rainbow.py."""
    import demo_rainbow
    assert True
