# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals, absolute_import
import numbers
import re

from ..pgraph import ParameterObject

re_split_port = re.compile(r'(.*?)/?\+')

class SimulationObject(ParameterObject):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #to help with debugging
        self.space_bonds = []
        self.autospaces = []

    def space_bond(self, *args, total_m = None):
        """
        TODO put this in the optics group, which doesn't have its own base class yet
        """
        from .. import optics
        last_object = None
        last_distance = None
        bond_list = []

        def dict_argument_add(arg):
            #doing some shitparsing here...
            nonlocal last_distance
            nonlocal bond_list
            nonlocal last_object

            item = arg['inj']
            item = item.strip()
            m = re_split_port.match(item)
            if not m:
                raise RuntimeError("'{}' doesn't look like a port reference".format(arg))
            objname = m.group(1)
            if '/' in objname:
                raise NotImplementedError("Cannot currently space_bond nested elements")
            if last_distance is None:
                raise RuntimeError("element injection (via dictionaries) only available if a distance is specified for previous element")
            else:
                if last_object is None:
                    raise RuntimeError("No previous object to connect space")
                space_name = "{}_{}/".format(last_object, objname)

                dist_fr = arg.get('fr', None)
                dist_to = arg.get('to', None)
                width   = arg.get('w', 0)
                if dist_to is not None:
                    if dist_fr is not None:
                        raise RuntimeError("cannot specify both from and to distances")
                    ins_distance = last_distance - dist_to - width
                    last_distance = dist_to
                else:
                    ins_distance = dist_fr
                    last_distance = last_distance - dist_fr - width

                Ssys = self[space_name] = optics.Space()
                Ssys['length[m]'] = ins_distance
                self.autospaces.append(space_name)

                #and inject the bond as well
                bond_list.append(space_name + '+A-t')
                bond_list.append(item)
            last_object = objname
            return

        def str_argument_add(arg):
            #doing some shitparsing here...
            nonlocal last_distance
            nonlocal bond_list
            nonlocal last_object
            spl = arg.split('|')
            if len(spl) > 1:
                for a in spl:
                    str_argument_add(a)
                return

            arg = arg.strip()
            m = re_split_port.match(arg)
            if not m:
                raise RuntimeError("'{}' doesn't look like a port reference".format(arg))
            objname = m.group(1)
            if '/' in objname:
                raise NotImplementedError("Cannot currently space_bond nested elements")
            if last_distance is None:
                bond_list.append(arg)
            else:
                if last_object is None:
                    raise RuntimeError("No previous object to connect space")
                space_name = "{}_{}/".format(last_object, objname)

                Ssys = self[space_name] = optics.Space()
                Ssys['length[m]'] = last_distance
                self.autospaces.append(space_name)

                #and inject the bond as well
                bond_list.append(space_name + '+A-t')
                bond_list.append(arg)
            last_object = objname
            last_distance = None
            return

        def argument_add(arg):
            nonlocal last_distance
            if isinstance(arg, str):
                str_argument_add(arg)
            elif isinstance(arg, dict):
                dict_argument_add(arg)
            elif isinstance(arg, numbers.Number):
                if last_distance is not None:
                    raise RuntimeError(
                        "Cannot have two sequential distances"
                    )
                last_distance = arg
            else:
                raise RuntimeError(
                    "Unrecognized space_bond type, must either be a"
                    " bond string, or a number representing a space"
                )
        for arg in args:
            argument_add(arg)
        self.bond_add(bond_list)
        self.space_bonds.append(bond_list)
        return bond_list


class OpticalObject(SimulationObject):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        with self._internal():
            self['label'] = 1

    def visit_mm_anno_label(self, pbg, view, descB):
        """
        """
        label = view['label']
        if isinstance(label, int):
            if label >= 0:
                return pbg.path_str_short(self, pre = label)
            else:
                return pbg.path_str(self, pre = label)
        else:
            return label

    def visit_mm_anno_description(self, pbg, view, descB):
        """
        """
        return None

