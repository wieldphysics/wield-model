# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import numpy as np
#from matplotlib.ticker import MultipleLocator, AutoMinorLocator

import collections
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
from matplotlib import legend_handler, lines
from matplotlib.text import OffsetFrom
import declarative

from ...optics import alm
from ...optics.alm.utils import (
    str_m,
)
from transient.utilities.mpl import mplfigB

from transient.utilities.mpl.stacked_plots import (
    generate_stacked_plot_ax,
)


class OverlapperPlotter(declarative.OverridableObject):
    fname     = None
    Z         = None
    N_points  = 300
    padding_m = None
    padding_rel = .00

    bbox_args = dict(
        alpha=0,
        ec = None,
        fc = None,
    )
    arrow_args = dict(
        arrowstyle="->",
        connectionstyle="angle,angleB=90,angleA=180,rad=3",
        linewidth = .5,
    )

    _overridable_object_save_kwargs = True
    _overridable_object_kwargs = None

    def __call__(self, *args, **kwargs):
        return self.plot(*args, **kwargs)

    def plot(
        self,
        overlapper      = None,
        fname           = None,
        Z               = None,
        use_in          = True,
        padding_m       = None,
        padding_rel     = None,
        ref0            = None,
        annotate = True,
        kw_target       = {},
        kwX_target      = {},
        kwY_target      = {},
        kwX             = dict(ls = '--', label = None),
        kwY             = dict(ls = '-'),
        kw              = dict(lw = 1),
        self_overlap    = False,
        reverse         = False,
        axB             = None,
        gouy_wrap       = 180,
        setup_xaxis     = True,
        setup_yaxis     = True,
        length_max_m    = None,
        length_max_in   = None,
        annotate_tags   = [],
        **kwargs
    ):
        fname = declarative.first_non_none(fname, self.fname)
        anno_target = overlapper.annotation_target()

        if use_in:
            z_unit_low = 1/.0254
            z_unit_top = 1
            if length_max_in is not None:
                length_max_m = length_max_in * .0254
        else:
            z_unit_low = 1
            z_unit_top = 1/.0254
            if length_max_m is None and length_max_in is not None:
                length_max_m = length_max_in * .0254

        Z           = declarative.first_non_none(Z, self.Z)
        padding_m   = declarative.first_non_none(padding_m, self.padding_m)
        padding_rel = declarative.first_non_none(padding_rel, self.padding_rel)

        if ref0 is not None:
            z0 = overlapper.object_z(ref0)
            print("Z0 is ", z0)
        else:
            z0 = 0

        if not reverse:
            z_reversed = None
            z_map = lambda z : z_unit_low * (z - z0)
        else:
            z_reversed = overlapper.length_m
            z_map = lambda z : (z_unit_low * (z_reversed - (z - z0)))

        if overlapper is None:
            raise RuntimeError("Must Provide a overlapper object to plot")

        if Z is None:
            if padding_m is None:
                padding_m = overlapper.length_m * padding_rel
            if length_max_m is not None:
                length_max_m = min(length_max_m, overlapper.length_m)
            else:
                length_max_m = overlapper.length_m
            Z = np.linspace(-padding_m, float(length_max_m) + padding_m, self.N_points)
            idx_pad_L = np.searchsorted(Z, 0)
            idx_pad_R = np.searchsorted(Z, float(overlapper.length_m))

        #axB = mplfigB(Nrows = 3)
        if axB is None:
            axB = generate_stacked_plot_ax(
                name_use_list = [
                    ('width', True),
                    ('iROC', True),
                    ('Gouy', True),
                    ('olap', self_overlap),
                ],
                width_phys_in = 10,
                heights_phys_in_default = 1.5,
                hspace = .08,
            )
        axB.width.set_ylabel(u'2$\\sigma$ intensity ⌀\n2w(z) [mm]')

        targetsB = dict()
        last_phase_X = 0
        last_phase_Y = 0
        last_zsub = 0
        #print(overlapper.target_list())
        for tname in overlapper.target_list():
            z_sub = Z
            #if not full_beams:
            #    if idx > 0:
            #        z_from = beam_zs[idx - 1]
            #        z_sub = z_sub[z_sub > z_from]
            #    if idx < len(beam_zs) - 1:
            #        z_to = beam_zs[idx + 1]
            #        z_sub = z_sub[z_sub < z_to]

            kw_line = dict(kw)
            kw_line.update(kw_target.get(tname, {}))
            kw_line.update(kwX)
            kw_line.update(kwX_target.get(tname, {}))

            zp_idx = np.searchsorted(z_sub, last_zsub)
            if zp_idx >= len(z_sub):
                zp_idx = -1
            last_zsub = z_sub[-1]

            qX = overlapper.z2target_qX(tname, z_sub)

            phase = np.unwrap(np.angle(qX.gouy_phasor))
            phase = phase - phase[zp_idx] + last_phase_X
            last_phase_X = phase[-1]

            Rline = axB.width.plot(z_map(z_sub), 2 * 1e3 * qX.W, **kw_line)
            line = Rline[0]
            color = line.get_color()
            kw_line.setdefault('color', color)
            axB.iROC.plot(z_map(z_sub), qX.R_inv, **kw_line)
            if gouy_wrap:
                axB.Gouy.plot(z_map(z_sub), (180 / np.pi * phase) % gouy_wrap, **kw_line)
            else:
                axB.Gouy.plot(z_map(z_sub), 180 / np.pi * phase, **kw_line)

            kw_line = dict(kw)
            if tname == anno_target:
                kw_line['label'] = "*" + tname
            else:
                kw_line['label'] = tname
            kw_line.update(kw_target.get(tname, {}))
            kw_line.update(kwY)
            kw_line.update(kwY_target.get(tname, {}))
            kw_line.setdefault('color', color)

            qY = overlapper.z2target_qY(tname, z_sub)

            phase = np.unwrap(np.angle(qY.gouy_phasor))
            phase = phase - phase[zp_idx] + last_phase_Y
            last_phase_Y = phase[-1]

            Rline = axB.width.plot(z_map(z_sub), 2 * 1e3 * qY.W, **kw_line)

            axB.iROC.plot(z_map(z_sub), qY.R_inv, **kw_line)
            if gouy_wrap:
                axB.Gouy.plot(z_map(z_sub), (180 / np.pi * phase) % gouy_wrap, **kw_line)
            else:
                axB.Gouy.plot(z_map(z_sub), 180 / np.pi * phase, **kw_line)

            if axB.olap:
                axB.olap.plot(z_map(z_sub), 100 * (1 - abs(qX.overlap_LG(qY))**2)/2, **kw_line)

        legend = axB.Gouy.legend(
            loc = 'upper left',
            ncol = 2,
            fontsize='medium',
        )
        if legend is not None:
            legend.get_frame().set_alpha(.9)

        axB.iROC.set_ylabel('iROC\n[1/m]')
        axB.Gouy.set_ylabel("Gouy Phase\n[deg]")
        if axB.olap:
            axB.olap.set_ylabel("Self X:Y Loss\nOverlap [% pwr]")

        if setup_xaxis:
            axB.ax_bottom.set_xlim(z_unit_low * min(Z - z0), z_unit_low * max(Z - z0))
            axB.ax_top_2 = axB.ax_top.twiny()
            axB.ax_top_2.set_xlim(z_unit_top * min(Z - z0), z_unit_top * max(Z - z0))

            if use_in:
                l = axB.ax_bottom.set_xlabel('Path [in]', labelpad = -8)
                l2 = axB.ax_top_2.set_xlabel('Path [m]', labelpad = -8)
            else:
                l = axB.ax_bottom.set_xlabel('Path [m]', labelpad = -8)
                l2 = axB.ax_top_2.set_xlabel('Path [in]', labelpad = -8)
            l.set_horizontalalignment('right')
            l.set_position((-0.01, +.1))
            l2.set_horizontalalignment('right')
            l2.set_position((-0.01, +.1))

        if annotate:
            descriptions = []
            descriptions += overlapper.plot_descriptions(
                reverse = reverse,
                tags = annotate_tags,
            )

            #cull the descriptions that extend past the plot
            descriptions2 = []
            for desc in descriptions:
                if desc.z1_m > length_max_m:
                    continue
                descriptions2.append(desc)
            descriptions = descriptions2
            #for tscname, tscfunc in transcribers.items():
            #    descriptions.update(tscfunc(overlapper.trans_center))
            self.annotate(
                axB = axB,
                descriptions = descriptions,
                use_in = use_in,
                z0 = z0,
                z_reversed = z_reversed,
                annotate = annotate,
                **kwargs
            )

            #for debugging
            #for descB in descriptions:
            #    axB.width.scatter(z_map(descB.z1_m), 1e3 * 2 * descB.q_start.W)
            #    axB.width.scatter(z_map(descB.z1_m), 1e3 * 2 * descB.q_start.W)

        xmin = z_unit_low * (0 - z0)
        xmax = z_unit_low * (float(overlapper.length_m) - z0)
        def limits_between(ax, xmin, xmax, spanscale_low = 1, spanscale_high = 1):
            lmax = -float('infinity')
            lmin = float('infinity')
            all_y = []
            for line in ax.get_lines():
                xd = line.get_xdata()
                yd = line.get_ydata()
                idx_pad_L = np.searchsorted(xd, xmin)
                idx_pad_R = np.searchsorted(xd, xmax)
                #only include if the line was in data coordinates
                if line.get_transform().contains_branch(ax.transData) and (idx_pad_R > idx_pad_L):
                    all_y.append(yd[idx_pad_L : idx_pad_R])
                    lmax = max(lmax, np.nanmax(yd[idx_pad_L : idx_pad_R]))
                    lmin = min(lmin, np.nanmin(yd[idx_pad_L : idx_pad_R]))
            ysorted = np.sort(np.concatenate(all_y))
            lmin = ysorted[int(len(ysorted) * .01)]
            lmax = ysorted[int(len(ysorted) * .99)]
            if not np.isfinite(lmax):
                lmax = None
            if not np.isfinite(lmin):
                lmin = None
            if lmin is not None and lmax is not None:
                _lmin = lmax - spanscale_high * (lmax - lmin)
                _lmax = lmin - spanscale_low * (lmin - lmax)
                lmin, lmax = _lmin, _lmax
            return lmin, lmax

        if setup_yaxis:
            axB.width.set_ylim(0, 1.1*limits_between(axB.width, xmin, xmax)[1])
            low, high = limits_between(axB.iROC, xmin, xmax, spanscale_low = 1.05, spanscale_high = 1.05)
            axB.iROC.set_ylim(low, high)
            #low, high = limits_between(axB.Gouy, xmin, xmax, spanscale_low = 1.05, spanscale_high = 1.05)
            #axB.Gouy.set_ylim(low, high)

        axB.finalize()
        axB.ax_bottom.minorticks_on()
        if setup_xaxis:
            axB.ax_top_2.minorticks_on()
        axB.width.minorticks_on()
        axB.iROC.minorticks_on()
        axB.Gouy.minorticks_on()

        axB.width.grid(which='minor', linewidth = 0.5, ls = ':')
        axB.iROC.grid(which='minor', linewidth = 0.5, ls = ':')
        axB.Gouy.grid(which='minor', linewidth = 0.5, ls = ':')

        axB.width.grid(which = 'major', linewidth = 1)
        axB.iROC.grid(which = 'major', linewidth = 1)
        axB.Gouy.grid(which = 'major', linewidth = 1)

        if fname is not None:
            axB.save(fname)
        return declarative.Bunch(locals())

    def annotate(
        self,
        axB,
        descriptions,
        use_in = False,
        z0 = 0,
        z_reversed = None,
        include_detuning = False,
        annotate = True,
    ):
        all_desc_by_z = []
        if use_in:
            z_unit_low = 1/.0254
        else:
            z_unit_low = 1

        if z_reversed is None:
            z_map = lambda z : z_unit_low * float(z)
        else:
            z_map = lambda z : z_unit_low * float(z_reversed - z)

        all_desc_by_z = list(descriptions)
        if not z_reversed:
            all_desc_by_z.sort(key = lambda d: (d['z_m'], d.get('L_m', float('inf'))))
        else:
            all_desc_by_z.sort(key = lambda d: (d['z_m'], -d.get('L_m', float('inf'))))

        zs = np.array([d['z_m'] for d in all_desc_by_z]) - z0

        xlow, xhigh = axB.ax_top.get_xlim()
        xmid = (xlow + xhigh)/2
        if annotate == 'top':
            idx_mid = len(zs)
        elif annotate == 'bottom':
            idx_mid = 0
        else:
            idx_mid = np.searchsorted(z_unit_low * (zs), xmid)

        if not z_reversed:
            left_list = all_desc_by_z[:idx_mid]
            right_list = all_desc_by_z[idx_mid:]
        else:
            left_list = all_desc_by_z[:idx_mid:-1]
            right_list = all_desc_by_z[idx_mid::-1]
        fsize_sep = 15

        if use_in:
            def desc_loc_format(z, L_m):
                if L_m > 0:
                    return u"{0:.3f}+{1:.2f}in".format(z * 100 / 2.54, L_m * 100 / 2.54)
                else:
                    return u"{0:.3f}in".format(z * 100 / 2.54)
        else:
            def desc_loc_format(z, L_m):
                if L_m > 0:
                    return u"{0}+{1}".format(str_m(z, 3, space=False), str_m(L_m, 2, space=False))
                else:
                    return u"{0}".format(str_m(z, 3, space=False))

        def desc_format(z, name, desc, L_m, left):
            locstr = desc_loc_format(z, L_m)
            if left:
                if name is not None:
                    desc = u"{0}:{name}: {desc}".format(locstr, desc = desc, name = name)
                else:
                    desc = u"{0}: {desc}".format(locstr, desc = desc)
            else:
                if name is not None:
                    desc = u"{desc}: {name}:{0}".format(locstr, desc = desc, name = name)
                else:
                    desc = u"{0}: {desc}".format(locstr, desc = desc)
            return desc

        def ready_desc(d):
            z_m     = d['z_m'] - z0
            width_m = d.get('L_m', None)
            span = d.get('span', True)
            line_kw = dict(
                lw = 0.5,
                ls = '--',
            )
            line_kw.update(d.get('line_kw', {}))
            anno_kw = dict(
                lw = 0.5,
            )
            anno_kw.update(d.get('anno_kw', {}))
            color   = d.get('color', None)

            if line_kw is None:
                line_kw = dict()

            if anno_kw is None:
                anno_kw = dict(ls = '--')
            if color is not None:
                anno_kw.setdefault('color', color)
                line_kw.setdefault('color', color)

            desc = desc_format(z_m, d.get('name'), d.get('desc'), d.get('L_m', 0), False)

            return declarative.Bunch(locals())

        for idx, d in enumerate(reversed(left_list)):
            d = ready_desc(d)
            #top elements
            arrowkw = dict(self.arrow_args)
            arrowkw.update(d.anno_kw)
            axB.ax_top.annotate(
                d.desc,
                xy=(z_map(d.z_m), 1), xycoords=axB.ax_top.get_xaxis_transform(),
                xytext=(0, 18 + fsize_sep*idx), textcoords=OffsetFrom(axB.ax_top.bbox, (1, 1), "points"),
                ha = "right", va = "bottom",
                bbox = self.bbox_args,
                arrowprops = arrowkw,
            )
            #axB.ax_top.annotate(
            #    desc,
            #    #'',
            #    xy=(0.0, 2), xytext=(0.0, 2),
            #    textcoords=OffsetFrom(an, (1, 1), "points"),
            #    ha="right", va="bottom",
            #    bbox=self.bbox_args,
            #)
            if d.span and d.width_m is not None and abs(d.width_m) > .001:
                axB.ax_top.annotate(
                    d.desc,
                    xy=(z_map(d.z_m + d.width_m), 1), xycoords=axB.ax_top.get_xaxis_transform(),
                    xytext=(0, 18 + fsize_sep*idx), textcoords=OffsetFrom(axB.ax_top.bbox, (1, 1), "points"),
                    ha = "right", va = "bottom",
                    bbox = self.bbox_args,
                    arrowprops = arrowkw,
                    alpha = 0,
                )
                axB.width.axvline(z_map(d.z_m), **d.line_kw)
                axB.width.axvline(z_map(d.z_m + d.width_m), **d.line_kw)
                axB.iROC.axvline(z_map(d.z_m), **d.line_kw)
                axB.iROC.axvline(z_map(d.z_m + d.width_m), **d.line_kw)
                axB.Gouy.axvline(z_map(d.z_m), **d.line_kw)
                axB.Gouy.axvline(z_map(d.z_m + d.width_m), **d.line_kw)
                if axB.olap:
                    axB.olap.axvline(z_map(d.z_m + d.width_m), **d.line_kw)
            else:
                axB.width.axvline(z_map(d.z_m), **d.line_kw)
                axB.iROC.axvline(z_map(d.z_m), **d.line_kw)
                axB.Gouy.axvline(z_map(d.z_m), **d.line_kw)
                if axB.olap:
                    axB.olap.axvline(z_map(d.z_m), **d.line_kw)

        for idx, d in enumerate(right_list):
            d = ready_desc(d)
            #bottom elements
            arrowkw = dict(self.arrow_args)
            arrowkw.update(d.anno_kw)
            axB.ax_top.annotate(
                d.desc,
                xy=(z_map(d.z_m), -.12), xycoords=axB.ax_bottom.get_xaxis_transform(),
                xytext=(0, -34 - fsize_sep*idx), textcoords=OffsetFrom(axB.ax_bottom.bbox, (0, 0), "points"),
                ha="left", va="bottom",
                bbox=self.bbox_args,
                arrowprops = arrowkw,
            )
            #axB.ax_top.annotate(
            #    '',#desc,
            #    xy=(0.0, 2),
            #    xytext=(0.0, 2),
            #    textcoords=OffsetFrom(an, (0, 0), "points"),
            #    ha="left", va="bottom",
            #    bbox=self.bbox_args,
            #)
            if d.span and d.width_m is not None and abs(d.width_m) > .001:
                axB.width.axvline(z_map(d.z_m), **d.line_kw)
                axB.width.axvline(z_map(d.z_m + d.width_m), **d.line_kw)
                axB.iROC.axvline(z_map(d.z_m), **d.line_kw)
                axB.iROC.axvline(z_map(d.z_m + d.width_m), **d.line_kw)
                axB.Gouy.axvline(z_map(d.z_m), **d.line_kw)
                axB.Gouy.axvline(z_map(d.z_m + d.width_m), **d.line_kw)
                if axB.olap:
                    axB.olap.axvline(z_map(d.z_m), **d.line_kw)
                    axB.olap.axvline(z_map(d.z_m + d.width_m), **d.line_kw)
                axB.ax_top.annotate(
                    d.desc,
                    xy=(z_map(d.z_m + d.width_m), -.12), xycoords=axB.ax_bottom.get_xaxis_transform(),
                    xytext=(0, -34 - fsize_sep*idx), textcoords=OffsetFrom(axB.ax_bottom.bbox, (0, 0), "points"),
                    ha="left", va="bottom",
                    bbox=self.bbox_args,
                    arrowprops = arrowkw,
                    alpha = 0,
                )
            else:
                axB.width.axvline(z_map(d.z_m), **d.line_kw)
                axB.iROC.axvline(z_map(d.z_m), **d.line_kw)
                axB.Gouy.axvline(z_map(d.z_m), **d.line_kw)
                if axB.olap:
                    axB.olap.axvline(z_map(d.z_m), **d.line_kw)
        return

    def plot_scan(
            self,
            overlapper,
            *args,
            axB = None,
            transverse = '',
            target_handles = None,
            black_handles = False,
            ncol = 4,
            **kwargs
    ):
        if transverse == '':
            axB_orig = axB
            if axB is None:
                axB = mplfigB(Ncols = 2, Nrows = 2, size_in = (12, 10))
                axB.legend_handles = []
                axB.legend_labels  = []
                if target_handles is None:
                    target_handles = collections.defaultdict(list)
                axB.target_handles = target_handles
            else:
                if target_handles is None:
                    target_handles = axB.target_handles

            #assign in axB to fill it with the points definitions
            P1B = self.plot_scan(
                overlapper = overlapper,
                axB = axB,
                axQ = axB.ax0,
                axLG = axB.ax2,
                transverse = 'x',
                target_handles = target_handles,
                **kwargs,
            )
            P2B = self.plot_scan(
                overlapper = overlapper,
                axQ = axB.ax1,
                axLG = axB.ax3,
                transverse = 'y',
                **kwargs,
            )
            ax = axB.ax1
            #box = ax.get_position()
            #ax.set_position([box.x0, box.y0 + box.height * 0.2,
            #                box.width, box.height * 0.8])

            # Put a legend below current axis
            if False and axB_orig is None:
                ax.legend(loc='upper left', bbox_to_anchor=(0.0, -0.12), ncol = 2)
            else:
                #TODO better legends
                handles = []
                labels = []
                longest = 1
                for (t, l), h in sorted(target_handles.items()):
                    if l is None:
                        continue
                    if black_handles:
                        if l == overlapper.target1:
                            h = mlines.Line2D([], [], ls='', marker='+', color = 'black'),
                        elif l == overlapper.target2:
                            h = mlines.Line2D([], [], ls='', marker='x', color = 'black'),
                        elif l == 'overlap {}\nwith {}'.format(overlapper.target1, overlapper.target2):
                            h = mlines.Line2D([], [], ls='', marker='o', color = 'black'),
                    h = [_ for _ in h if _ is not None]
                    if len(h) > longest:
                        longest = len(h)
                    if len(h) > 1:
                        handles.append(tuple(h))
                    else:
                        handles.append(h[0])
                    labels.append(str(l))
                    print(labels[-1], handles[-1])
                ax.legend(
                    handles, labels,
                    loc='upper left',
                    bbox_to_anchor=(-0.05, -0.12),
                    ncol = ncol,
                    handler_map={
                        tuple: HandlerTuple(ndivide=None),
                        lines.Line2D: HandlerLine2Dv(),
                    },
                    handlelength = longest * .8,

                )
            #axB.ax1.legend(framealpha = 1)
            return axB
        else:
            return self._plot_scan(
                overlapper,
                *args,
                axB = axB,
                transverse = transverse,
                target_handles = target_handles,
                **kwargs
            )
    def _plot_scan(
        self,
        overlapper = None,
        axB = None,
        axQ = None,
        axLG = None,
        fname = None,
        include_contour = True,
        include_limits = True,
        scans = None,
        color = None,
        label1 = None,
        label2 = None,
        labelT1 = None,
        labelT2 = None,
        transverse = 'x',
        s1_split = None,
        s2_split = None,
        no_lines = False,
        target_handles = None,
        group_full = False,
        labelgroup = None,
    ):
        transverse = transverse.lower()
        assert(transverse in ['x', 'y', ''])

        label1_orig  = label1
        label2_orig  = label2

        if labelT1 is not None:
            labelT1 = overlapper.target1
        labelT1_orig = labelT1
        if labelT2 is not None:
            labelT2 = overlapper.target2
        labelT2_orig = labelT2

        if axQ is None and axLG is None:
            if axB is None:
                axB = mplfigB(Ncols = 2, size_in = (12, 5))
            axQ = axB.ax0
            axLG = axB.ax1
            axQ.grid(b=True, alpha = .2)
            axLG.grid(b=True, alpha = .2)

        calc = overlapper.compile_overlap_calculation()
        pbg_orig = overlapper.pbg
        pbg = pbg_orig.copy()

        for scanD in scans:
            scanD = declarative.Bunch(scanD)
            num = scanD.get('num', 1)
            val_abs = scanD.get('values', None)
            val_rel = scanD.get('values_rel', None)
            val_shift = scanD.get('values_shift', None)

            oval = pbg_orig.get_parameter(scanD.parameter)

            if val_abs is None:
                vals = np.asarray(oval)
            else:
                vals = np.asarray(val_abs)

            if val_shift is not None:
                vals = vals + np.asarray(val_shift)

            if val_rel is not None:
                vals = vals + oval * val_rel

            if num == 1:
                vals_M = vals.reshape(-1, 1)
                if s1_split is None:
                    s1_split = len(vals.reshape(-1)) // 2
            elif num == 2:
                vals_M = vals.reshape(1, -1)
                if s2_split is None:
                    s2_split = len(vals.reshape(-1)) // 2
            else:
                raise RuntimeError("Can only handle up to 2 scans")
            pbg.override_value(scanD.parameter, vals_M)

        if s1_split is None:
            s1_split = 0
        if s2_split is None:
            s2_split = 0

        qB = calc.calculate_Qs(pbg)

        def plot_points(
                ax,
                X,
                Y,
                color = None,
                label = 'test',
                marker = 'o',
                s = 16,
        ):
            nonlocal label1
            nonlocal label2
            points = None
            line1 = None
            line2 = None
            if X.shape == () or (X.shape[0] == 1 and X.shape[1] == 1):
                path = ax.scatter(
                    X,
                    Y,
                    s = s,
                    marker = marker,
                    color = color,
                    label = label,
                )
                if color is None:
                    color = path.get_edgecolor()[0, :3]
                points = path
                line1 = None
            elif X.shape[0] == 1:
                path = ax.scatter(
                    X,
                    Y,
                    s = s,
                    marker = marker,
                    color = color,
                    label = label,
                )
                points = path
                if color is None:
                    color = path.get_edgecolor()[0, :3]
                if not no_lines:
                    path, = ax.plot(
                        X[0, s2_split:],
                        Y[0, s2_split:],
                        color = color,
                        lw = 1,
                        ls = '--',
                        label = label1,
                    )
                    #clear label to prevent repeats in legend
                    label1 = None
                    line1 = path
                    path, = ax.plot(
                        X[0, :s2_split+1],
                        Y[0, :s2_split+1],
                        color = color,
                        ls = '--',
                        lw = 0.5,
                    )
            elif X.shape[1] == 1:
                path = ax.scatter(
                    X,
                    Y,
                    s = s,
                    marker = marker,
                    color = color,
                    label = label,
                )
                points = path
                if color is None:
                    color = path.get_edgecolor()[0, :3]
                if not no_lines:
                    path, = ax.plot(
                        X[s1_split:],
                        Y[s1_split:],
                        color = color,
                        ls = '-',
                        lw = 1,
                        label = label1,
                    )
                    #clear label to prevent repeats in legend
                    label1 = None
                    line1 = path
                    path, = ax.plot(
                        X[:s1_split+1],
                        Y[:s1_split+1],
                        ls = '-',
                        color = color,
                        lw = 0.5,
                    )
            else:
                path = ax.scatter(
                    X,
                    Y,
                    s = s,
                    marker = marker,
                    color = color,
                    label = label,
                )
                points = path
                if color is None:
                    color = path.get_edgecolor()[0, :3]
                if not no_lines:
                    path, = ax.plot(
                        X[s1_split:, s2_split],
                        Y[s1_split:, s2_split],
                        color = color,
                        lw = 1,
                        ls = '-',
                        label = label1,
                    )
                    label1 = None
                    line1 = path
                    path, = ax.plot(
                        X[:s1_split+1, s2_split],
                        Y[:s1_split+1, s2_split],
                        color = color,
                        ls = '-',
                        lw = 0.5,
                    )
                    path, = ax.plot(
                        X[s1_split, s2_split:],
                        Y[s1_split, s2_split:],
                        color = color,
                        lw = 1,
                        ls = '--',
                        label = label2,
                    )
                    label2 = None
                    line2 = path
                    path, = ax.plot(
                        X[s1_split, :s2_split + 1],
                        Y[s1_split, :s2_split + 1],
                        color = color,
                        ls = '--',
                        lw = .5,
                    )
            return declarative.Bunch(
                color = color,
                line1 = line1,
                line2 = line2,
                points = points,
            )

        if axQ is not None:
            if transverse == 'x':
                if overlapper.target1 is not None:
                    pointsB1 = plot_points(
                        ax = axQ,
                        X = 1/qB.t1qX.R,
                        Y = qB.t1qX.W * 1e3,
                        color = color,
                        label = labelT1,
                        s = 25,
                        marker = '+',
                    )
                    labelT1 = None
                    color = pointsB1.color
                if overlapper.target2 is not None:
                    pointsB2 = plot_points(
                        ax = axQ,
                        X = 1/qB.t2qX.R,
                        Y = qB.t2qX.W * 1e3,
                        color = color,
                        label = labelT2,
                        marker = 'x',
                    )
                    labelT2 = None
                    color = pointsB2.color
            else:
                if overlapper.target1 is not None:
                    pointsB1 = plot_points(
                        ax = axQ,
                        X = 1/qB.t1qY.R,
                        Y = qB.t1qY.W * 1e3,
                        color = color,
                        label = labelT1,
                        s = 25,
                        marker = '+',
                    )
                    labelT1 = None
                    color = pointsB1.color
                if overlapper.target2 is not None:
                    pointsB2 = plot_points(
                        ax = axQ,
                        X = 1/qB.t2qY.R,
                        Y = qB.t2qY.W * 1e3,
                        color = color,
                        label = labelT2,
                        marker = 'x',
                    )
                    labelT2 = None
                    color = pointsB2.color

        if axLG is not None:
            #olap_00X, olap_02X = qB.t1qX.overlap_LG_2mode(qB.t2qX)
            #olap_00Y, olap_02Y = qB.t1qY.overlap_LG_2mode(qB.t2qY)
            #olap_02 = olap_00Y*olap_02X + olap_00X*olap_02Y
            #axLG.scatter(olap_02.real, olap_02.imag)

            if overlapper.target1 is not None and overlapper.target2:
                if transverse == 'x':
                    olap_00, olap_02 = qB.t1qX.overlap_LG_2mode(qB.t2qX)
                else:
                    olap_00, olap_02 = qB.t1qY.overlap_LG_2mode(qB.t2qY)
                pointsB12 = plot_points(
                    ax = axLG,
                    X = olap_02.real,
                    Y = olap_02.imag,
                    color = color,
                    s = 16,
                    marker = '.',
                )
                color = pointsB12.color

        if include_limits:
            if axQ is not None:
                mmB_Qx = data_minmax(axQ)
                axQ.set_xlim(mmB_Qx.xmin, mmB_Qx.xmax)
                axQ.set_ylim(mmB_Qx.ymin, mmB_Qx.ymax)

            if axLG is not None:
                mmB_LG = data_minmax(axLG)
                if mmB_LG.rmax is not None:
                    axLG.set_xlim(-mmB_LG.rmax, mmB_LG.rmax)
                    axLG.set_ylim(-mmB_LG.rmax, mmB_LG.rmax)
                axLG.set_aspect(1)

        if include_contour:
            if axQ is not None:
                Qref = None
                if transverse == 'x':
                    if overlapper.target2 is not None and qB.t2qX.W.shape == ():
                        Qref = qB.t2qX
                    elif overlapper.target1 is not None and qB.t1qX.W.shape == ():
                        Qref = qB.t1qX
                else:
                    if overlapper.target2 is not None and qB.t2qY.W.shape == ():
                        Qref = qB.t2qY
                    elif overlapper.target1 is not None and qB.t1qY.W.shape == ():
                        Qref = qB.t1qY

                if Qref is not None:
                    xmin, xmax = axQ.get_xlim()
                    ymin, ymax = axQ.get_ylim()
                    X_iR, Y_W = np.meshgrid(
                        np.linspace(xmin, xmax, 200),
                        np.linspace(ymin, ymax, 200) / 1e3,
                    )
                    Qplot = alm.ComplexBeamParam.from_W_R(Y_W, 1/X_iR, wavelength_m = 1064e-9)
                    Lolap = 1 - abs(Qref.overlap_LG(Qplot))**2
                    CS = axQ.contour(
                        X_iR, Y_W * 1e3, 100 * Lolap,
                        levels = [.1, .3, 1, 3, 10, 30],
                        colors = ['silver', 'gray', 'black', 'purple', 'red', 'red'],
                        alpha = .7,
                        linewidths = 1,
                    )
                    axQ.clabel(CS, inline=1, fontsize=8)

            if axLG is not None:
                xmin, xmax = axLG.get_xlim()
                ymin, ymax = axLG.get_ylim()
                X, Y = np.meshgrid(
                    np.linspace(xmin, xmax, 200),
                    np.linspace(ymin, ymax, 200),
                )
                Lolap = X**2 + Y**2
                CS = axLG.contour(
                    X,
                    Y,
                    100 * Lolap,
                    levels = [.1, .3, 1, 3, 10, 30],
                    colors = ['silver', 'gray', 'black', 'purple', 'red', 'red'],
                    alpha = .7,
                    linewidths = 1,
                )
                axLG.clabel(CS, inline=1, fontsize=8)

        if axQ is not None:
            axQ.set_ylabel("Beam Radius w [mm]")
            axQ.set_xlabel("inverse ROC [Diopters]")
            if transverse == 'x':
                axQ.set_title("Beam Parameter Plot [X]")
            else:
                axQ.set_title("Beam Parameter Plot [Y]")
            axQ.grid(b=True, alpha = .2)
        if axLG is not None:
            axLG.set_ylabel("LG1 imaginary")
            axLG.set_xlabel("LG1 real")
            if transverse == 'x':
                axLG.set_title("Rel. mismatch LG0->LG1 (X parameter coupling)")
            else:
                axLG.set_title("Rel. mismatch LG0->LG1 (Y parameter coupling)")
            axLG.grid(b=True, alpha = .2)

        if fname is not None:
            axB.save(fname)
        handlesB = declarative.Bunch()
        handlesB.B1  = pointsB1
        handlesB.B2  = pointsB2
        handlesB.B12 = pointsB12
        if axB is None:
            axB = declarative.Bunch()
        axB.handlesB = handlesB

        if target_handles is not None:
            if label1_orig is not None:
                if handlesB.B1.line1 is not None:
                    target_handles[4, label1_orig].append(
                        handlesB.B1.line1
                    )
                elif handlesB.B2.line1 is not None:
                    target_handles[4, label1_orig].append(
                        handlesB.B2.line1
                    )
                elif handlesB.B12.line1 is not None:
                    target_handles[4, label1_orig].append(
                        handlesB.B12.line1
                    )
            if label2_orig is not None:
                if handlesB.B1.line2 is not None:
                    target_handles[5, label2_orig].append(
                        handlesB.B1.line2
                    )
                elif handlesB.B2.line2 is not None:
                    target_handles[5, label2_orig].append(
                        handlesB.B2.line2
                    )
                elif handlesB.B12.line2 is not None:
                    target_handles[5, label2_orig].append(
                        handlesB.B12.line2
                    )
            if handlesB.B1.points is not None:
                target_handles[1, overlapper.target1].append(
                    handlesB.B1.points
                )
                if group_full:
                    target_handles[4, label1_orig].append(
                        handlesB.B1.points
                    )
                    target_handles[5, label2_orig].append(
                        handlesB.B1.points
                    )
            if handlesB.B2.points is not None:
                target_handles[2, overlapper.target2].append(
                    handlesB.B2.points
                )
                if group_full:
                    target_handles[4, label1_orig].append(
                        handlesB.B2.points
                    )
                    target_handles[5, label2_orig].append(
                        handlesB.B2.points
                    )
            if handlesB.B12.points is not None:
                target_handles[3, 'overlap {}\nwith {}'.format(overlapper.target1, overlapper.target2)].append(
                    handlesB.B12.points
                )
                if group_full:
                    target_handles[4, label1_orig].append(
                        handlesB.B12.points
                    )
                    target_handles[5, label2_orig].append(
                        handlesB.B12.points
                    )
            if labelgroup is not None:
                grp = target_handles[6, labelgroup]
                grp.extend(
                    h for h in [
                        handlesB.B1.line1,
                        handlesB.B1.line2,
                        handlesB.B1.points,
                        handlesB.B2.points,
                        handlesB.B12.points,
                    ] if h is not None
                )
        return axB


def data_minmax(
    ax,
    spanscale = 1.2,
    Xspanscale = None,
    Yspanscale = None,
    Xspanscale_high = None,
    Xspanscale_low  = None,
    Yspanscale_high = None,
    Yspanscale_low  = None,
):
    if Xspanscale is None:
        Xspanscale = spanscale
    if Xspanscale_low is None:
        Xspanscale_low = Xspanscale
    if Xspanscale_high is None:
        Xspanscale_high = Xspanscale
    if Yspanscale is None:
        Yspanscale = spanscale
    if Yspanscale_low is None:
        Yspanscale_low = Yspanscale
    if Yspanscale_high is None:
        Yspanscale_high = Yspanscale

    ymax = -float('infinity')
    ymin = float('infinity')
    all_y = []
    xmax = -float('infinity')
    xmin = float('infinity')
    all_x = []
    rmax = 0
    all_r = []
    for line in ax.get_children():
        import matplotlib as mpl
        if isinstance(line, mpl.collections.Collection):
            try:
                off = line.get_offsets()
                xd, yd = off.T
                trans = line.get_offset_transform()
            except AttributeError:
                continue
        elif isinstance(line, mpl.lines.Line2D):
            try:
                xd = line.get_xdata()
                yd = line.get_ydata()
                trans = line.get_transform()
            except AttributeError:
                continue
        #only include if the line was in data coordinates
        if trans.contains_branch(ax.transData):
            all_y.append(yd)
            ymax = max(ymax, np.nanmax(yd))
            ymin = min(ymin, np.nanmin(yd))
            all_x.append(xd)
            xmax = max(xmax, np.nanmax(xd))
            xmin = min(xmin, np.nanmin(xd))
            rd = (xd**2 + yd**2)**0.5
            all_r.append(rd)
            rmax = max(rmax, np.nanmax(rd))
    ysorted = np.sort(np.concatenate(all_y))
    ymin = ysorted[int(len(ysorted) * .01)]
    ymax = ysorted[int(len(ysorted) * .99)]
    xsorted = np.sort(np.concatenate(all_x))
    xmin = xsorted[int(len(xsorted) * .01)]
    xmax = xsorted[int(len(xsorted) * .99)]
    rsorted = np.sort(np.concatenate(all_r))
    rmax = rsorted[int(len(rsorted) * .99)]
    if not np.isfinite(ymax):
        ymax = None
    if not np.isfinite(ymin):
        ymin = None
    if not np.isfinite(xmax):
        xmax = None
    if not np.isfinite(xmin):
        xmin = None
    if not np.isfinite(rmax):
        rmax = None
    if ymin is not None and ymax is not None:
        _ymin = ymax - Yspanscale_high * (ymax - ymin)
        _ymax = ymin - Yspanscale_low * (ymin - ymax)
        ymin, ymax = _ymin, _ymax
    if xmin is not None and xmax is not None:
        _xmin = xmax - Xspanscale_high * (xmax - xmin)
        _xmax = xmin - Xspanscale_low * (xmin - xmax)
        xmin, xmax = _xmin, _xmax
    if rmax is not None:
        rmax = spanscale * rmax
    return declarative.Bunch(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        rmax = rmax,
    )


class HandlerLine2Dv(legend_handler.HandlerLine2D):
    """
    Handler for `.Line2D` instances.
    """

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize,
                       trans):

        xdata, xdata_marker = self.get_xdata(legend, xdescent, ydescent,
                                             width, height, fontsize)
        xdata = np.asarray(xdata)

        ydata = np.linspace(-0.5 * height, 1.5 * height, xdata.shape[0])
        legline = lines.Line2D(xdata, ydata)

        self.update_prop(legline, orig_handle, legend)
        legline.set_drawstyle('default')
        legline.set_marker("")

        legline_marker = lines.Line2D(xdata_marker, ydata[:len(xdata_marker)])
        self.update_prop(legline_marker, orig_handle, legend)
        legline_marker.set_linestyle('None')
        if legend.markerscale != 1:
            newsz = legline_marker.get_markersize() * legend.markerscale
            legline_marker.set_markersize(newsz)
        # we don't want to add this to the return list because
        # the texts and handles are assumed to be in one-to-one
        # correspondence.
        legline._legmarker = legline_marker

        legline.set_transform(trans)
        legline_marker.set_transform(trans)

        return [legline, legline_marker]

