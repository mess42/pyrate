#!/usr/bin/python3
import optical_system
import ray
import numpy as np
import matplotlib.pyplot as plt


# TODO: Method headers, Comments
# TODO: three spot diagrams: vertex spot diagram, chief ray spot diagram, centroid spot diagram


def drawLayout2d(ax, s, list_of_raypaths, list_of_raypath_colors = 0):
    if list_of_raypath_colors == 0:
        list_of_raypath_colors = map(lambda e: 'blue', range(len(list_of_raypaths)))
    map(lambda r: r[0].draw2d(s, ax, color=r[1]), zip(list_of_raypaths, list_of_raypath_colors))
    s.draw2d(ax)


def drawSpotDiagram(ax, s, raypath, surfaceno):
    # perhaps: raypath should give points at certain surface number back
    # such that we respect privacy of the raypath
    #print raypath.raybundles[surfaceno].o
    pts = raypath.raybundles[surfaceno].o
    ax.plot(pts[0], pts[1],'.')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')

    # xlabel, ylabel auf spot beziehen

    rmsspot = raypath.raybundles[surfaceno].getRMSspotSizeCentroid()

    ax.text(0.05, 0.05,'Centroid RMS spot radius: '+str(1000.*rmsspot)+' um', transform=ax.transAxes)

    if surfaceno < 0:
        surfaceno = len(s.surfaces) - 1
    ax.set_title('Spot diagram at surface no. '+str(surfaceno))
