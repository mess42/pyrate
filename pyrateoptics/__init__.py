#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2018
               by     Moritz Esslinger moritz.esslinger@web.de
               and    Johannes Hartung j.hartung@gmx.net
               and    Uwe Lippmann  uwe.lippmann@web.de
               and    Thomas Heinze t.heinze@uni-jena.de
               and    others

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

# Use this file for convenience functions which should be called
# at the main package level.

import logging
from distutils.version import StrictVersion

import numpy as np

import matplotlib
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider

from .core.iterators import OptimizableVariableKeyIterator
from .raytracer.analysis.optical_system_analysis import OpticalSystemAnalysis
from .raytracer.optical_system import OpticalSystem
from .raytracer.localcoordinates import LocalCoordinates
from .raytracer.optical_element import OpticalElement
from .raytracer.surface import Surface
from .raytracer.surface_shape import accessible_shapes, LinearCombination
from .raytracer.globalconstants import (numerical_tolerance,
                                        standard_wavelength,
                                        degree)
from .raytracer.ray import RayBundle, RayPath
from .raytracer.material.material_isotropic import ConstantIndexGlass
from .raytracer.material.material_glasscat import GlassCatalog

# TODO: provide convenience classes for building a builduplist which could be
# transferred to the build...functions


def build_rotationally_symmetric_optical_system(builduplist, **kwargs):
    """
    Convenience function to build up a centrosymmetric system with
    conic lenses.

    :param builduplist: (list of tuple)
             elements are (r, cc, thickness, mat, name, optdict)
             r - radius of curvature (float)
             cc - conic constant (float)
             thickness - (float)
             mat - material index or name (str)
                   if convertable to float,
                       a ConstantIndexGlass will be created
                   if not, a material from
                       the refractiveindex.info will be created
            name - name of surf (str)
            optdict - {"is_mirror": True|False, "is_stop": True|False}

    :param material_db_path: (str)
             path to the refractiveindex.info yml-database

    :return (s, stdseq): (tuple)
             s is an OpticalSystem object
             stdseq is a sequence for sequential raytracing
    """

    builduplist_build_simple_os = []
    for (r, cc, thickness, mat, name, optdict) in builduplist:
        curv = 0
        if abs(r) > numerical_tolerance:
            curv = 1./r
        else:
            curv = 0.
        coordbrkdict = {"decz": thickness}
        surfdict = {"shape": "Conic", "curv": curv, "cc": cc}
        builduplist_build_simple_os.append((surfdict, coordbrkdict, mat, name,
                                            optdict))

    return build_simple_optical_system(builduplist_build_simple_os, **kwargs)


def build_simple_optical_element(lc0, builduplist, material_db_path="",
                                 name=""):
    """
    Convenience function to build up an optical element.

    """
    logger = logging.getLogger(__name__)

    elem = OpticalElement.p(lc0, name=name)

    refname = lc0.name
    lastmat = None
    surflist_for_sequence = []

    gcat = GlassCatalog(material_db_path)

    for (surfdict, coordbreakdict, mat, surf_name, optdict) in builduplist:
        lc = elem.addLocalCoordinateSystem(
            LocalCoordinates.p(name=surf_name + "_lc", **coordbreakdict),
            refname=refname)
        shapetype = "shape_" + surfdict.pop("shape", "Conic")
        aperture = surfdict.pop("aperture", None)
        if shapetype == "shape_LinearCombination":
            # linear combination accepts pairs of coefficients and surfdicts

            list_of_coefficients_and_shapes =\
                surfdict.get("list_of_coefficients_and_shapes", [])
            new_list_coeffs_shapes = []
            for (ind, (coeff_part, surfdict_part)) in\
                    enumerate(list_of_coefficients_and_shapes, 1):
                shapetype_part = "shape_" + surfdict_part.pop("shape", "Conic")
                new_list_coeffs_shapes.append((coeff_part,
                                               accessible_shapes[
                                                   shapetype_part].p
                                               (lc,
                                                name=name + "_shape" +
                                                str(ind),
                                                **surfdict_part)))

            actsurf = Surface.p(lc, name=surf_name + "_surf",
                                aperture=aperture,
                                shape=LinearCombination.p(lc,
                                                          name=surf_name +
                                                          "_linearcombi",
                                                          list_of_coefficients_and_shapes=new_list_coeffs_shapes))
        else:
            actsurf = Surface.p(lc, name=surf_name + "_surf",
                                aperture=aperture,
                                shape=accessible_shapes[shapetype].p
                                (lc, name=name + "_shape", **surfdict))
        logger.debug("mat=%s" % repr(mat))
        # TODO: evaluation function which adds a material depending on type
        # of argument and can handle:
        # * a string (database),
        # * an object (direct material definition),
        # * a floating (later maybe complex number) (constant index glass)
        if mat is not None:
            use_floating_point_value_for_constant_index_glass = False
            try:
                n = float(mat)
            except ValueError:
                use_floating_point_value_for_constant_index_glass = False
            else:
                use_floating_point_value_for_constant_index_glass = True
            if isinstance(mat, str) and not use_floating_point_value_for_constant_index_glass:
                gcat.get_material_dictFromLongName(mat)
                elem.addMaterial(mat,
                                 gcat.createGlassObjectFromLongName(lc, mat))
            elif use_floating_point_value_for_constant_index_glass:
                mat = "constantindexglass_" + str(mat)
                elem.addMaterial(mat, ConstantIndexGlass.p(lc, n=n))

        elem.addSurface(surf_name, actsurf, (lastmat, mat))
        logger.info("Added surface: %s at material boundary %s" % (surf_name,
                                                                   (lastmat,
                                                                    mat)))

        lastmat = mat
        refname = lc.name
        surflist_for_sequence.append((surf_name, optdict))

    return (elem, (name, surflist_for_sequence))


def build_simple_optical_system(builduplist, material_db_path="", name=""):

    """
    Convenience function to build up system with simple lenses.

    :param builduplist: (list of tuple of dicts)
            elements are (surfdict, coordbreakdict, mat, name, optdict)
            surfdict - {"shape": "Conic, ...", "curv": ..., "cc": ...}
            (such that surfdict can be used for **kwargs in surface_shape)
            coordbreakdict - {"decz": thickness, decx: ..., decy: ...,
                              tiltx: ..., ..., order: 0 or 1}
            (such that coordbreakdict can be used for **kwargs in
             LocalCoordinates)
            mat - material index or name (str)
            name - name of surf (str)
            optdict - {"is_mirror":True|False, "is_stop":True|False}

    :param material_db_path: (str)
             path to the refractiveindex.info yml-database

    :return (s, stdseq): (tuple)
             s is an OpticalSystem object
             stdseq is a sequence for sequential raytracing
    """
    logger = logging.getLogger(__name__)
    logger.info("Creating simple optical system")
    s = OpticalSystem.p(name=name)

    lc0 = s.addLocalCoordinateSystem(
        LocalCoordinates.p(name="object", decz=0.0),
        refname=s.rootcoordinatesystem.name)

    elem_name = "stdelem"
    logger.info("Element name %s" % (elem_name,))

    (elem, elem_seq) = build_simple_optical_element(
        lc0, builduplist,
        material_db_path=material_db_path, name=elem_name)
    s.addElement(elem_name, elem)

    s.material_background.set_name("background")
    stdseq = [(elem_seq)]
    logger.info("Created simple optical system")

    return (s, stdseq)


def build_optical_system(builduplist, material_db_path="", name=""):
    """
    Convenience function to build up non-centrosymmetric optical system.

    TODO: parameters and return type
    """
    logger = logging.getLogger(__name__)
    logger.info("Creating multiple element optical system")
    s = OpticalSystem(name=name)
    lc0 = s.addLocalCoordinateSystem(
        LocalCoordinates(name="object", decz=0.0),
        refname=s.rootcoordinatesystem.name)

    full_elements_seq = []
    refname = lc0.name
    for (element_list, coordbreakdict, elem_name) in builduplist:
        logger.info("Element name %s" % (elem_name,))
        lc = s.addLocalCoordinateSystem(
            LocalCoordinates(name=elem_name + "_lc", **coordbreakdict),
            refname=refname)
        refname = lc.name

        (elem, elem_seq) = build_simple_optical_element(
            lc, element_list,  # builduplist?
            material_db_path=material_db_path, name=elem_name)
        full_elements_seq.append(elem_seq)
        s.addElement(elem_name, elem)

    s.material_background.set_name("background")
    logger.info("Created multiple element optical system")

    return (s, full_elements_seq)


def draw(os, rays=None,
         hold_on=False,
         do_not_draw_surfaces=[],
         do_not_draw_raybundles=[],
         interactive=False,
         show_box=True,
         linewidth=1.0,
         export_type="pdf",
         export=None, **kwargs):

    """
    Convenience function for drawing optical system and list of raybundles
    Use figsize=(..,..) to control aspect ratio for export.

    :param os - OpticalSystem
    :param rb - list of raybundles

    """

    if export is not None:
        interactive = False

    axis_color = "lightgoldenrodyellow"

    dpi = kwargs.pop("dpi", None)
    figsize = kwargs.pop("figsize", None)

    fig = plt.figure(1, dpi=dpi, figsize=figsize)

    if not show_box:
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
    else:
        ax = fig.add_subplot(1, 1, 1)

    if interactive:
        fig.subplots_adjust(left=0.25, bottom=0.25)

        xz_angle_slider_size = [0.25, 0.15, 0.65, 0.03]
        up_angle_slider_size = [0.25, 0.10, 0.65, 0.03]

        if StrictVersion(matplotlib.__version__) < StrictVersion("2.0.0"):
            xz_angle_slider_ax = fig.add_axes(xz_angle_slider_size,
                                              axisbg=axis_color)
            up_angle_slider_ax = fig.add_axes(up_angle_slider_size,
                                              axisbg=axis_color)
        else:
            xz_angle_slider_ax = fig.add_axes(xz_angle_slider_size,
                                              facecolor=axis_color)
            up_angle_slider_ax = fig.add_axes(up_angle_slider_size,
                                              facecolor=axis_color)

        xz_angle_slider = Slider(xz_angle_slider_ax,
                                 "XZ angle",
                                 0.0, 360.0,
                                 valinit=0.0, valfmt="%0.0f")
        up_angle_slider = Slider(up_angle_slider_ax,
                                 "UP angle",
                                 0.0, 360.0,
                                 valinit=0.0, valfmt="%0.0f")

    ax.axis("equal")
    if StrictVersion(matplotlib.__version__) < StrictVersion("2.0.0"):
        ax.set_axis_bgcolor("white")
    else:
        ax.set_facecolor("white")

    def draw_rays(ax, rays, do_not_draw_raybundles=[], **kwargs):
        if rays is not None:
            if isinstance(rays, list):
                for rpl in rays:
                    ray_color = tuple(np.random.random(3))
                    if isinstance(rpl, list):
                        for rp in rpl:
                            ray_color = tuple(np.random.random(3))
                            rp.draw2d(ax, color=ray_color,
                                      do_not_draw_raybundles=do_not_draw_raybundles,
                                      **kwargs)
                    elif isinstance(rpl, tuple):
                        (rl, ray_color) = rpl
                        if isinstance(rl, list):
                            # draw(s, [([rp1, ..], color1), (....)])
                            for r in rl:
                                r.draw2d(ax, color=ray_color, **kwargs)
                        else:
                            # draw(s, [(rp1, color1), (....)])
                            rl.draw2d(ax, color=ray_color, **kwargs)
                    else:
                        rpl.draw2d(ax, color=ray_color,
                                   do_not_draw_raybundles=do_not_draw_raybundles,
                                   **kwargs)
            elif isinstance(rays, RayPath):
                # draw(s, raypath)
                ray_color = tuple(np.random.random(3))
                rays.draw2d(ax, color=ray_color,
                            do_not_draw_raybundles=do_not_draw_raybundles,
                            **kwargs)
            elif isinstance(rays, RayBundle):
                # draw(s, raybundle)
                ray_color = tuple(np.random.random(3))
                rays.draw2d(ax, color=ray_color, **kwargs)
            elif isinstance(rays, tuple):
                (rl, ray_color) = rays
                if isinstance(rl, list):
                    # draw(s, ([raypath1, ...], color))
                    for r in rl:
                        r.draw2d(ax, color=ray_color, **kwargs)
                else:
                    # draw(s, (raypath, color))
                    rl.draw2d(ax, color=ray_color,
                              do_not_draw_raybundles=do_not_draw_raybundles,
                              **kwargs)

    def sliders_on_changed(val):
        ax.clear()

        round_val_xz = np.round(xz_angle_slider.val)
        round_val_up = np.round(up_angle_slider.val)

        phi = round_val_xz*degree
        theta = round_val_up*degree

#        new_ex = np.array([np.cos(phi),
#                           0,
#                           np.sin(phi)])
#        new_up = np.array([-np.sin(theta),
#                           np.cos(theta),
#                           0.])
        new_ex = np.array([np.cos(phi)*np.cos(theta),
                           np.sin(theta),
                           np.sin(phi)*np.cos(theta)])
        new_up = np.array([-np.cos(phi)*np.sin(theta),
                           np.cos(theta),
                           -np.sin(phi)*np.sin(theta)])
        inyzplane = np.abs(round_val_xz) < numerical_tolerance\
            and np.abs(round_val_up) < numerical_tolerance
        os.draw2d(ax, color="grey", inyzplane=inyzplane,
                  plane_normal=new_ex, up=new_up, **kwargs)
        draw_rays(ax, rays, plane_normal=new_ex, up=new_up, **kwargs)

    if interactive:
        xz_angle_slider.on_changed(sliders_on_changed)
        up_angle_slider.on_changed(sliders_on_changed)

    draw_rays(ax, rays, linewidth=linewidth,
              do_not_draw_raybundles=do_not_draw_raybundles, **kwargs)
    os.draw2d(ax, color="grey", linewidth=linewidth,
              do_not_draw_surfaces=do_not_draw_surfaces, **kwargs)

    if export is not None:
        # ax.autoscale(enable=True, axis='both', tight=True)
        fig.savefig(export, format=export_type,
                    bbox_inches='tight', pad_inches=0)

    if not hold_on:
        plt.show()


def raytrace(s, seq, numrays, rays_dict, bundletype="collimated",
             traceoptions={}, wave=standard_wavelength):
    """
    Convenience function for raytracing.
    """
    osa = OpticalSystemAnalysis(s, seq)
    osa.aim(numrays, rays_dict, bundletype=bundletype, wave=wave)

    return osa.trace(**traceoptions)


def listOptimizableVariables(os, filter_status=None, max_line_width=None):
    """
    Convenience function to list all optimizable variables within
    an object.

    :param os (ClassWithOptimizableVariables) -
            Container with optimizable variables.
    :param filter_status (str or None) -
            filter table for "variable", "fixed", ...

    :returns os.getAllVariables()
    """

    lst = OptimizableVariableKeyIterator(os).variables_dictionary

    def shorten_string(s, maxlen=None, intermediate_string="..."):

        if maxlen is None:
            return s
        else:
            if len(s) > maxlen:
                to_remove = len(s) - maxlen + len(intermediate_string)
                interpos = (len(s) - to_remove) // 2
                return (s[:interpos] + intermediate_string +
                        s[len(s) - interpos:])
            else:
                return s

    def print_table(table):
        col_width = [max(len(str(x)) for x in col) for col in zip(*table)]
        for line in table:
            print(" ".join("{:{}}".format(x, col_width[i])
                           for i, x in enumerate(line)))

    table = [(shorten_string(a, maxlen=max_line_width),
              b.var_type(),
              str(b.evaluate()))
             for (a, b) in sorted(lst.items(), key=lambda x: (
                 len(x[0].split('.')) - 1, x[0]))]
    # sort by number of . in dict key and afterwards by lexical order

    if filter_status is not None:
        table = [(a, b, c) for (a, b, c) in table if b == filter_status]

    print_table(table)

    return lst
