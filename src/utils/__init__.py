"""Shared helpers hoisted out of the figure notebooks.

Layout follows the scanpy convention used by Moorman et al., Nature 2024:
`pl` plotting, `tl` tools.

Only helpers with a single unambiguous definition live here. Helpers that exist in
several diverging versions (`makeneighborgraph` and `joinsets`, which differ between
`SAND.ipynb` and `tl/runSAND_transcript.py`, and `p_to_star`) are deliberately left
inline in the notebooks and cells that use them, so that each panel keeps the version
it was actually run with.
"""
