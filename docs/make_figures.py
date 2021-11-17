import os.path
from typing import Tuple

import numpy as np
import matplotlib.pyplot as pp
import mpl_toolkits.mplot3d.art3d as art3d

from xopto import pf

DOCS_DIR = os.path.dirname(os.path.abspath(__file__))
EPS = np.finfo(np.float).eps

def radial_accumulator(ax=None, start: float = 0.4, stop=1.0, n: int = 6,
                       logscale: bool = False):
    start, stop, n = float(start), float(stop), int(n)
    dr = (stop - start)/n

    if logscale:
        rs = np.exp(np.linspace(np.log(start), np.log(stop), n + 1))
    else:
        rs = np.linspace(start, stop, n + 1)

    for r2, r1 in ((rs[-1], rs[-2]), (rs[1], rs[0])):
        c1 = pp.Circle((0.0, 0.0), r2, fill=True, color='lightgray')
        c2 = pp.Circle((0.0, 0.0), r1, fill=True, color='white')
        ax.add_artist(c1)
        ax.add_artist(c2)

    fig = None
    if ax is None:
        fig, ax = pp.subplots()

    for r in rs:
        if r > 0.0:
            circ = pp.Circle((0.0, 0.0), r,
                             clip_on=False, color='k', fill=False, ls='-')
            ax.add_artist(circ)

    px = py = 0.05
    if start > EPS:
        pp.plot([-px, px], [0, 0], '-k')
        pp.plot([0, 0], [-py, py], '-k')
        ax.text(0, -px,
                '$position = (x, y)$', ha='center', va='top')

        xy = -start*np.cos(np.pi/4), start*np.sin(np.pi/4)
        pp.annotate('', xy=xy, xytext=(0, 0), arrowprops=dict(arrowstyle="-|>"))
        ax.text(px*2 + xy[0], xy[1], 'start', ha='left', va='center')

    ax.text(px + stop/2**0.5, -(py + stop/2**0.5),
           '$n$', ha='left', va='top')

    xy = stop*np.cos(np.pi/4), stop*np.sin(np.pi/4)
    pp.annotate('', xy=xy, xytext=(0, 0), arrowprops=dict(arrowstyle="-|>"))
    ax.text(xy[0] + px, xy[1] + py, 'stop')

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_aspect('equal', 'box')
    ax.set_axis_off()

    return fig, ax


def cartesian_accumulator(ax=None, start=-0.8, stop=0.8, nx=10, ny=10):
    start, stop, nx, ny = float(start), float(stop), int(nx), int(ny)

    dx, dy = (stop - start)/nx, (stop - start)/ny

    fig = None
    if ax is None:
        fig, ax = pp.subplots()

    for x, y in ((start, start), (stop - dx, start), (stop - dx, stop - dy)):
        rect = pp.Rectangle((x, y), dx, dy, color='lightgray', fill=True)
        ax.add_artist(rect)

    for x in np.linspace(start, stop, nx + 1):
        pp.plot([x, x], [start, stop], '-k')
    for y in np.linspace(start, stop, ny + 1):
        pp.plot([start, stop], [y,y], '-k')

    py = -0.05
    pp.annotate('', xy=(stop, start + py), xytext=(start, start + py),
                arrowprops=dict(arrowstyle="-|>"))
    pp.text(start, start + 2*py, '$start_x$', va='top', ha='left')
    pp.text(stop, start + 2*py, '$stop_x$', va='top', ha='right')
    pp.text(0, start + 2*py, '$n_x$', va='top', ha='center')

    px = 0.05
    pp.annotate('', xy=(stop + px, start), xytext=(stop + px, stop),
                arrowprops=dict(arrowstyle="-|>"))
    pp.text(stop + 2*px, stop, '$start_y$', va='top',   ha='left')
    pp.text(stop + 2*px, start, '$stop_y$', va='bottom',ha='left')
    pp.text(stop + 2*px, 0.0, '$n_y$', va='center', ha='left')

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_aspect('equal', 'box')
    ax.set_axis_off()

    return fig, ax

def symmetricx_accumulator(ax=None, width : float = 0.8,
                           n_half: int = 10, logscale: bool = False):

    fig = None
    if ax is None:
        fig, ax = pp.subplots()

    width, n = float(width), int(n_half)
    height = 0.8

    if logscale:
        x2 = np.exp(np.linspace(np.log(EPS), np.log(width), n_half + 1))
        x1 = -x2[::-1]
    else:
        x1 = np.linspace(-width, 0, n_half + 1)
        x2 = np.linspace(0.0, width, n_half + 1)

    rect = pp.Rectangle((x1[0], -height), x1[1] - x1[0], 2*height, color='lightgray', fill=True)
    ax.add_artist(rect)
    rect = pp.Rectangle((x1[-2], -height), x1[-1] - x1[-2], 2*height, color='lightgray', fill=True)
    ax.add_artist(rect)
    rect = pp.Rectangle((x2[0], -height), x2[1] - x2[0], 2*height, color='lightgray', fill=True)
    ax.add_artist(rect)
    rect = pp.Rectangle((x2[-2], -height), x2[-1] - x2[-2], 2*height, color='lightgray', fill=True)
    ax.add_artist(rect)

    for x in x1:
        pp.plot([x, x], [-0.8, 0.8], '-k')
    for x in x2:
        pp.plot([x, x], [-0.8, 0.8], '-k')

    py = -0.05
    pp.annotate('', xy=(width, -height + py), xytext=(0.0, -height + py),
                arrowprops=dict(arrowstyle="-|>"))
    pp.annotate('', xy=(-width, -height + py), xytext=(0.0, -height + py),
                arrowprops=dict(arrowstyle="-|>"))
    pp.text(x1[0], -height + 2*py, '$center - range$', va='top', ha='left')
    pp.text(-width*0.5, -height - 0.1 + 2*py, '$n_{half}$', va='top', ha='center')
    pp.text(width, -height + 2*py, '$center + range$', va='top', ha='right')
    pp.text(width*0.5, -height - 0.1 + 2*py, '$n_{half}$', va='top', ha='center')
    pp.text(0, -height + 2*py, '$center$', va='top', ha='center')

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_aspect('equal', 'box')
    ax.set_axis_off()

    return fig, ax

def fluence_rz_accumulator(ax=None, start=-0.8, stop=0.8, nr=10, nz=10):
    start, stop, nx, ny = float(start), float(stop), int(nr), int(nz)

    dx, dy = (stop - start)/nx, (stop - start)/ny

    fig = None
    if ax is None:
        fig, ax = pp.subplots()

    for x, y in ((start, start), (stop - dx, start), (start, stop - dy)):
        rect = pp.Rectangle((x, y), dx, dy, color='lightgray', fill=True)
        ax.add_artist(rect)

    for x in np.linspace(start, stop, nx + 1):
        pp.plot([x, x], [start, stop], '-k')
    for y in np.linspace(start, stop, ny + 1):
        pp.plot([start, stop], [y,y], '-k')

    py = -0.05
    pp.annotate('', xy=(stop, start + py), xytext=(start, start + py),
                arrowprops=dict(arrowstyle="-|>"))
    pp.text(start, start + 2*py, '$start_r$', va='top', ha='left')
    pp.text(stop, start + 2*py, '$stop_r$', va='top', ha='right')
    pp.text(0, start + 2*py, '$n_r$', va='top', ha='center')

    px = 0.05
    pp.annotate('', xy=(start - px, start), xytext=(start - px, stop),
                arrowprops=dict(arrowstyle="-|>"))
    pp.text(start - 2*px, stop,  '$start_z$', va='top',    ha='right')
    pp.text(start - 2*px, start, '$stop_z$',  va='bottom', ha='right')
    pp.text(start - 2*px, 0.0,   '$n_z$',     va='center', ha='right')

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_aspect('equal', 'box')
    ax.set_axis_off()

    return fig, ax

def fluence_rz_accumulator_3d(ax=None, zstart=0.0, zstop=2.0,
                              rstart=0.0, rstop=1.0, nz=10, nr=10):
    zstart, zstop, nr = float(zstart), float(zstop), int(nz)
    rstart, rstop, nr = float(rstart), float(rstop), int(nr)

    dz, dr = (zstop - zstart)/nr, (rstop - rstart)/nz

    fig = None
    if ax is None:
        fig, ax = pp.subplots()

    zs = np.linspace(zstart, zstop, nz + 1)
    rs = np.linspace(rstart, rstop, nr + 1)

    for z in zs[1:-1]:
        circ = pp.Circle((0.0, 0.0), rs[-1], ls=':', fill=False, edgecolor='black', facecolor='white')
        ax.add_patch(circ)
        art3d.pathpatch_2d_to_3d(circ, z=z, zdir="z")

    circ = pp.Circle((0.0, 0.0), rs[-1], fill=False, edgecolor='black', facecolor='white')
    ax.add_patch(circ)
    art3d.pathpatch_2d_to_3d(circ, z=zstart, zdir="z")

    #circ = pp.Circle((0.0, 0.0), rs[-1], fill=True, facecolor='white', zorder=0)
    #ax.add_patch(circ)
    #art3d.pathpatch_2d_to_3d(circ, z=zstop, zdir="z")
    for r in rs[::-1]:
        circ = pp.Circle((0.0, 0.0), r, fill=False, edgecolor='black', facecolor='white')
        ax.add_patch(circ)
        art3d.pathpatch_2d_to_3d(circ, z=zstop, zdir="z")


    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(0, 2)
    #ax.set_aspect('equal', 'box')
    ax.set_axis_off()

    return fig, ax


def six_around_one(ax, dcore: float, dcladding: float, diameter: float,
                   spacing: float = None, xyrange: Tuple[float, float] = None):
    fig = None
    if ax is None:
        fig, ax = pp.subplots()

    if spacing is None:
        spacing = dcladding

    # probe
    circ = pp.Circle((0.0, 0.0), diameter*0.5, color='gray', fill=True)
    ax.add_artist(circ)

    # central source fiber
    circ = pp.Circle((0.0, 0.0), dcladding*0.5, color='blue', fill=True)
    ax.add_artist(circ)
    circ = pp.Circle((0.0, 0.0), dcore*0.5, color='yellow', fill=True)
    ax.add_artist(circ)

    # surrounding fibers
    for index in np.arange(6):
        center = spacing*np.cos(index*np.pi/3), spacing*np.sin(index*np.pi/3)
        circ = pp.Circle(center, dcladding*0.5, color='blue', fill=True)
        ax.add_artist(circ)
        circ = pp.Circle(center, dcore*0.5, color='lightblue', fill=True)
        ax.add_artist(circ)

    if xyrange is not None:
        ax.set_xlim(*xyrange)
        ax.set_ylim(*xyrange)

    ax.set_aspect('equal', 'box')
    ax.set_axis_off()

    return fig, ax

def linear_array(ax, n:int, dcore: float, dcladding: float, diameter: float,
                 spacing: float = None, xyrange: Tuple[float, float] = None):
    fig = None
    if ax is None:
        fig = pp.figure()
        ax = fig.add_subplot(projection='3d')

    if spacing is None:
        spacing = dcladding

    # probe
    circ = pp.Circle((0.0, 0.0), diameter*0.5, color='gray', fill=True)
    ax.add_artist(circ)

    # fibers - from left to right
    for index in np.arange(n):
        core_color = 'yellow' if n > 1 and index == 0 else 'lightblue'
        center = -0.5*(n - 1)*spacing + index*spacing, 0.0
        circ = pp.Circle(center, dcladding*0.5, color='blue', fill=True)
        ax.add_artist(circ)
        circ = pp.Circle(center, dcore*0.5, color=core_color, fill=True)
        ax.add_artist(circ)

    if xyrange is not None:
        ax.set_xlim(*xyrange)
        ax.set_ylim(*xyrange)

    ax.set_aspect('equal', 'box')
    ax.set_axis_off()

    return fig, ax


def pf_hg(ax, g: tuple or list):
    fig = None
    if ax is None:
        fig, ax = pp.subplots()

    theta = np.linspace(0, np.pi, 1000)
    theta_deg = np.rad2deg(theta)
    ct = np.cos(theta)
    for g_ in g:
        pp.semilogy(theta_deg, pf.Hg(g_)(ct), label='$g=${:.1f}'.format(g_))

    ax.set_xlabel('Scattering angle $\Theta$ (°)')
    #ax.set_title('Henyey-Greenstein (HG) phase function')
    pp.legend()

    return fig, ax

def pf_mhg(ax, g: tuple or list, b: tuple or list):
    fig = None
    if ax is None:
        fig, ax = pp.subplots()

    theta = np.linspace(0, np.pi, 1000)
    theta_deg = np.rad2deg(theta)
    ct = np.cos(theta)
    for g_ in g:
        for b_ in b:
            pp.semilogy(theta_deg, pf.MHg(g_, b_)(ct),
                        label='$g=${:3.1f}, $\\beta=${:3.1f}'.format(g_, b_))

    ax.set_xlabel('Scattering angle $\Theta$ (°)')
    #ax.set_title('Modified Henyey-Greenstein (MHG) phase function')
    pp.legend()

    return fig, ax

def pf_gk(ax, g: tuple or list, a: tuple or list):
    fig = None
    if ax is None:
        fig, ax = pp.subplots()

    theta = np.linspace(0, np.pi, 1000)
    theta_deg = np.rad2deg(theta)
    ct = np.cos(theta)
    for g_ in g:
        for a_ in a:
            pp.semilogy(theta_deg, pf.Gk(g_, a_)(ct),
                        label='$g=${:3.1f}, $\\alpha=${:3.1f}'.format(g_, a_))

    ax.set_xlabel('Scattering angle $\Theta$ (°)')
    #ax.set_title('Modified Henyey-Greenstein (MHG) phase function')
    pp.legend()

    return fig, ax

def pf_mie(ax, ns: float, nm: float, d: tuple or list, w: float):
    fig = None
    if ax is None:
        fig, ax = pp.subplots()

    theta = np.linspace(0, np.pi, 1000)
    theta_deg = np.rad2deg(theta)
    ct = np.cos(theta)
    for d_ in d:
        pp.semilogy(theta_deg, pf.Mie(ns, nm, d_, w)(ct),
                    label='$diameter=${:.2f} μm'.format(d_*1e6))

    ax.set_xlabel('Scattering angle $\Theta$ (°)')
    #ax.set_title('Mie phase function of spherical particles')
    pp.legend()

    return fig, ax

def export(fig, filename):
    fig.savefig(os.path.join(DOCS_DIR, 'source', 'static', filename + '.svg'))


if __name__ == '__main__':
    SHOW = True

    fig, ax = pp.subplots()
    pf_hg(ax, g=(0.1, 0.5, 0.9))
    export(fig, 'mcml_comparison_hg')

    fig, ax = pp.subplots()
    pf_hg(ax, g=(0.1, 0.3, 0.5, 0.7, 0.9))
    export(fig, 'mcml_hg')

    fig, ax = pp.subplots()
    pf_mhg(ax, g=(0.1, 0.5, 0.9), b=(0.0, 0.4, 1.0))
    export(fig, 'mcml_mhg')

    fig, ax = pp.subplots()
    pf_gk(ax, g=(0.1, 0.5, 0.9), a=(-0.5, 0.5, 1.5))
    export(fig, 'mcml_gk')

    fig, ax = pp.subplots()
    pf_mie(ax, 1.462, 1.337, [0.25e-6, 0.5e-6, 1.0e-6, 2.0e-6, 4.0e-6], 500e-9)
    export(fig, 'mcml_mie_fused_silica')

    fig, ax = pp.subplots()
    pf_mie(ax, 1.603, 1.337, [0.25e-6, 0.5e-6, 1.0e-6, 2.0e-6, 4.0e-6], 500e-9)
    export(fig, 'mcml_mie_polystyrene')

    fig = pp.figure(figsize=(6.4*1.5, 4.8))
    ax = fig.add_subplot(1, 2, 1)
    fluence_rz_accumulator(ax)

    #fig = pp.figure()
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    fluence_rz_accumulator_3d(ax)
    export(fig, 'mcfluence_fluence_rz_3d')

    fig, ax = pp.subplots()
    radial_accumulator(ax)
    export(fig, 'mcdetector_radial')

    fig, ax = pp.subplots()
    radial_accumulator(ax, start=0)
    export(fig, 'mcdetector_radial_start_0')

    fig, ax = pp.subplots()
    radial_accumulator(ax, logscale=True)
    export(fig, 'mcdetector_radial_logscale')

    fig, ax = pp.subplots()
    radial_accumulator(ax, start=EPS, logscale=True, n=100)
    export(fig, 'mcdetector_radial_logscale_start_0')

    fig, ax = pp.subplots()
    cartesian_accumulator(ax)
    export(fig, 'mcdetector_cartesian')

    fig, ax = pp.subplots()
    symmetricx_accumulator(ax)
    export(fig, 'mcdetector_symmetricx')

    fig, ax = pp.subplots()
    symmetricx_accumulator(ax, logscale=1, n_half=100)
    export(fig, 'mcdetector_symmetricx_logscale')

    for dcore, dcladding in ((200e-6, 220e-6), (400e-6, 420e-6)):
        fig, ax = pp.subplots()
        six_around_one(ax, dcore, dcladding, 6e-3, xyrange=[-3.5e-3, 3.5e-3])
        export(fig, 'mcsurface_six-around-one-{:.0f}um'.format(dcore*1e6))

    fig, ax = pp.subplots(1, 2)
    for ax, dcore, dcladding in ((ax[0], 200e-6, 220e-6), (ax[1], 400e-6, 420e-6)):
        six_around_one(ax, dcore, dcladding, 6e-3, xyrange=[-3.5e-3, 3.5e-3])
    export(fig, 'mcsurface_six-around-one-200-400um')

    for dcore, dcladding, n, filename in (
            (200e-6, 220e-6, 6, 'mcsurface_six-linear-200um'),
            (100e-6, 120e-6, 1, 'mcsurface_single-fiber-100um'),
            (200e-6, 220e-6, 1, 'mcsurface_single-fiber-200um'),
            (400e-6, 420e-6, 1, 'mcsurface_single-fiber-400um'),
            (800e-6, 820e-6, 1, 'mcsurface_single-fiber-800um')):
        fig, ax = pp.subplots()
        linear_array(ax, n, dcore, dcladding, 6e-3, xyrange=[-3.5e-3, 3.5e-3])
        export(fig, filename)

    fig, ax = pp.subplots(2, 2)
    for ax, dcore, dcladding in (
            (ax[0, 0], 100e-6, 120e-6),
            (ax[0, 1], 200e-6, 220e-6),
            (ax[1, 0], 400e-6, 420e-6),
            (ax[1, 1], 800e-6, 820e-6)):
        linear_array(ax, 1, dcore, dcladding, 6e-3, xyrange=[-3.0e-3, 3.0e-3])
    export(fig, 'mcsurface_single-fiber-100-200-400-800um')

    if SHOW:
        pp.show()

    exit(0)