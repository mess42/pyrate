"""
Pyrate - Optical raytracing based on Python
"""

def test_smoke_doublet():
    """Smoke test based on demo_doublet.py."""

    # we use the matplotlib's do nothing backend for testing
    # matplotlib/lib/matplotlib/backends/backend_template.py
    import matplotlib
    matplotlib.use('Template')

    import demo_doublet
    assert False 
