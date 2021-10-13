# -*- coding: utf-8 -*-
"""
"""
from wavestate.bunch import Bunch
from ... import optics
from ...optics.alm.utils import (
    str_m,
)


def transporter2objects(algo_mm, q, transporter):
    """
    Creates a list of "descriptions", which are Bunches containing a set of data
    for annotations
    """
    descriptions = []
    mat2 = mat1 = transporter.inc_build_mat[0]
    idx_last = 0
    z_last = 0

    def object_desc(obj, descB):
        view = algo_mm.pbg.view(obj)
        label = obj.visit_mm_anno_label(algo_mm.pbg, view, descB)
        desc = obj.visit_mm_anno_description(algo_mm.pbg, view, descB)
        if desc is None:
            desc = label
        else:
            desc = '{} {}'.format(label, desc)
        descB.desc = desc
        descB.object = obj
    for idx, (ol_fr, ol_to) in enumerate(path_pairs(transporter.oLp_path)):
        idx_fr = transporter.inc_ol2idx.get(ol_fr, None)
        idx_to = transporter.inc_ol2idx.get(ol_to, None)

        #most of the transitions should occur from an element to itself, but at the edges,
        #we can add the start and end objects which occur at a transition between objects using this special case
        if ol_fr[0] != ol_to[0]:
            #print(ol_fr, ol_to, idx_fr, idx_to)
            if idx == 0:
                descB = Bunch(
                    z_m = 0,
                    L_m = 0,
                    z1_m = 0,
                    z2_m = 0,
                    color = 'purple',
                    q_start = q,
                    q_end   = q,
                    mat1 = mat1,
                    mat2 = mat2,
                )
                object_desc(ol_fr[0], descB)
                descriptions.append(descB)
            elif idx == len(transporter.oLp_path) - 2:
                if idx_fr is None:
                    continue
                z_m = transporter.inc_build_len[idx_fr]
                mat1 = transporter.inc_build_mat[idx_fr]
                q_end = q.propagate_matrix(mat1)
                descB = Bunch(
                    z_m = z_m,
                    L_m = 0,
                    z1_m = z_m,
                    z2_m = z_m,
                    color = 'purple',
                    q_start = q_end,
                    q_end   = q_end,
                    mat1 = mat1,
                    mat2 = mat1,
                )
                object_desc(ol_to[0], descB)
                descriptions.append(descB)
            #ignore links between objects
            continue

        if idx_fr is None or idx_to is None:
            z1 = z2 = z_last
            mat2 = mat1 = transporter.inc_build_mat[idx_last]
        else:
            if idx_to - idx_fr <= 0:
                continue

            z1 = transporter.inc_build_len[idx_fr]
            z2 = transporter.inc_build_len[idx_to]
            z_last = z2
            idx_last = idx_to

            mat1 = transporter.inc_build_mat[idx_fr]
            mat2 = transporter.inc_build_mat[idx_to]

        descB = Bunch(
            z_m = z1,
            z1_m = z1,
            z2_m = z2,
            L_m = z2 - z1,
            color = 'purple',
            q_start = q.propagate_matrix(mat1),
            q_end   = q.propagate_matrix(mat2),
            mat1 = mat1,
            mat2 = mat2,
        )
        object_desc(ol_fr[0], descB)
        descriptions.append(descB)
    return descriptions


def annotate_tags(pbg, descriptions, tags):
    for tag in tags:
        obj = tag['obj']
        if isinstance(obj, str):
            obj = pbg.resolve_object(tag['obj'])
        for idx, desc in enumerate(descriptions):
            if desc['object'] == obj:
                break
        else:
            raise RuntimeError("Tag object not found")
        to = tag.get('to', None)
        fr = tag.get('fr', None)
        if to is not None:
            desc2 = descriptions[idx - 1]
            assert(desc2)
            assert(isinstance(desc2.object, optics.Space))
            assert(desc2.L_m > to)
            q = desc.q_start.propagate_distance(-to)
            insert = slice(idx - 1, idx - 1)
            z_m = desc.z1_m - to
        elif fr is not None:
            desc2 = descriptions[idx + 1]
            assert(desc2)
            assert(isinstance(desc2.object, optics.Space))
            assert(desc2.L_m > fr)
            q = desc.q_end.propagate_distance(fr)
            insert = slice(idx, idx)
            z_m = desc.z2_m + fr
        else:
            raise RuntimeError("Tag needs either to or from component")

        label = tag['label'].format(
            **q.str_kw()
        )

        descriptions[insert] = [Bunch(
            z_m     = z_m,
            z1_m    = z_m,
            z2_m    = z_m,
            L_m     = 0,
            span    = False,
            object  = None,
            desc    = label,
            color   = tag.get('color', 'red'),
            q_start = q,
            q_end   = q,
        )]
    return descriptions

def annotate_waists(descriptions):
    descriptions2 = []
    for desc in descriptions:
        if isinstance(desc.object, optics.Space):
            #desc.color = 'green'
            #p = self.mm.pbg.view(desc.object)
            if -desc.q_start.Z < desc.L_m and -desc.q_start.Z >= 0:
                q_start = desc.q_start.propagate_distance(-desc.q_start.Z)
                z_m = desc.z_m - desc.q_start.Z
                descriptions2.append(Bunch(
                    z_m = z_m,
                    z1_m = z_m,
                    z2_m = z_m,
                    L_m = 0,
                    desc = 'waist Ø:{}, ZR:{}'.format(str_m(2 * q_start.W0), str_m(q_start.ZR)),
                    color = 'green',
                    object = None,
                    q_start = q_start,
                    q_end   = q_start,
                ))
            descriptions2.append(desc)
        else:
            descriptions2.append(desc)
    return descriptions2


def annotate_qspaces(descriptions, reverse = False):
    """
    Takes a description list, locates the spaces
    """
    descriptions2 = []
    idx = 0
    while idx < len(descriptions):
        desc = descriptions[idx]
        descriptions2.append(desc)
        idx += 1
        if isinstance(desc.object, optics.Space):
            if not desc.object.annotate_as_space:
                continue
            q_start = desc.q_start
            L_m = desc.L_m
            z_m = desc.z_m
            z1_m = desc.z1_m
            z2_m = desc.z2_m
            if L_m == 0:
                continue
            if reverse:
                q_start = q_start.propagate_distance(L_m)
                z_m = z_m + L_m
            descriptions2.append(Bunch(
                z_m     = z_m,
                z1_m    = z1_m,
                z2_m    = z2_m,
                L_m     = L_m,
                span    = False,
                object  = None,
                desc    = 'space q:{}'.format(q_start.str_short()),
                color   = 'orange',
                q_start = desc.q_start,
                q_end   = desc.q_end,
            ))
    return descriptions2


def annotate_clearspaces(descriptions):
    descriptions2 = []
    for desc in descriptions:
        if isinstance(desc.object, optics.Space):
            if not desc.object.annotate_as_space:
                descriptions2.append(desc)
        else:
            descriptions2.append(desc)
    return descriptions2


def path_pairs(link_path):
    last_link = link_path[0]
    for next_link in link_path[1:]:
        yield last_link, next_link
        last_link = next_link
    return
