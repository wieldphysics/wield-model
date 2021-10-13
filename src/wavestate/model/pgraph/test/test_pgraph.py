# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals, absolute_import

from wavestate.model import pgraph

def test_obj_build():
    obj_sys = pgraph.ParameterObject()
    obj_sys['obj1/'] = pgraph.ParameterObject()
    obj_sys['obj2/'] = pgraph.ParameterObject()
    obj_sys['obj1/imed/obj_i/'] = pgraph.ParameterObject()

    #obj_sys['obj1/imed/'] = pgraph.ParameterObject()
    obj1 = obj_sys['obj1/']
    obj1['imed/'] = pgraph.ParameterObject()

    print(obj_sys._reference_dict)

    pg = pgraph.ParameterGraph(obj_sys)
    print(pg.object_path_built)

    print('----------')
    print(pg.object_parameters)
    pg.print_parameters()


def test_obj_parameters():
    obj_sys = pgraph.ParameterObject()
    obj_sys['obj1/'] = pgraph.ParameterObject()
    obj_sys['obj2/'] = pgraph.ParameterObject()
    obj_sys['obj1/obj3/'] = pgraph.ParameterObject()

    obj_sys['param1'] = 100
    obj_sys['param2'] = 200

    obj_sys.set_assign('param1A', 'param1', lambda x : x + 1)
    obj_sys.set_assign('param1B', 'obj1/paramA', lambda x : x + 1)
    obj_sys.set_assign('obj1/paramB', 'obj1/paramA', lambda x : x + 1)
    obj_sys.set_assign('obj2/paramC', 'param1A', lambda x : x + 1)

    @obj_sys.deco_parameter
    def param_many(p):
        return p['obj2/paramC'] + p['obj1/paramB'] + p['param1A']

    @obj_sys.deco_multi_parameter(assignments = ['X1', 'X2', 'obj2/X3'])
    def param_many(p):
        p['X1'] = p['obj2/paramC'] + p['obj1/paramB'] + p['param1A']
        p['X2'] = 'test'
        p['obj2/X3'] = 'test'

    obj1 = obj_sys['obj1/']
    obj1['imed/'] = pgraph.ParameterObject()
    obj1['paramA'] = 300
    obj1['paramB'] = -300

    print(obj_sys._reference_dict)

    pg = pgraph.ParameterGraph(obj_sys)
    print(pg.object_path_built)

    print('----------')
    pg.print_parameters_eval()
    assert(pg.dict_parameters_eval() == {
        '/param1': 100,
        '/param_many': 504,
        '/X1': 504,
        '/param2': 200,
        'obj2/paramC': 102,
        'obj2/X3': 'test',
        '/param1A': 101,
        '/param1B': 301,
        'obj1/paramB': 301,
        'obj1/paramA': 300,
        '/X2': 'test'
    })


def test_obj_coverage():
    obj_sys = pgraph.ParameterObject()
    obj_sys['obj1/'] = pgraph.ParameterObject()
    obj_sys['obj2/'] = pgraph.ParameterObject()
    obj_sys['obj1/obj3/'] = pgraph.ParameterObject()
    obj_sys.set_assign('obj2/obj3/', 'obj1/obj3/')

    obj_sys['p1'] = 100

    obj1 = obj_sys['obj1/']
    obj1['imed/'] = pgraph.ParameterObject()
    obj1['obj3/paramA'] = 300
    obj2 = obj_sys['obj2/']

    #obj2['obj3/paramB'] = -300
    obj_sys['obj2/obj3/paramB'] = -300

    print(obj_sys._reference_dict)

    pg = pgraph.ParameterGraph(obj_sys)
    print(pg.object_path_built)

    print('----------')
    pg.print_parameters_eval()
    #print(pg.dict_parameters_eval())
    assert(pg.dict_parameters_eval() == {
        '/p1': 100,
        'obj1/obj3/paramB': -300,
        'obj1/obj3/paramA': 300,
    })


def test_obj_coverage2():
    sys = pgraph.ParameterObject()
    sys['o/'] = obj_sys = pgraph.ParameterObject()
    obj_sys['obj1/'] = pgraph.ParameterObject()
    obj_sys['obj2/'] = pgraph.ParameterObject()
    obj_sys['obj1/obj3/'] = pgraph.ParameterObject()
    #obj_sys.set_assign('obj3/', 'obj1/obj3/')
    obj_sys['obj3/'] = 'obj1/obj3/'

    obj_sys['p1'] = 100

    obj1 = obj_sys['obj1/']
    obj1['imed/'] = pgraph.ParameterObject()
    obj1['obj3/paramA'] = 300
    obj2 = obj_sys['obj2/']

    #obj2['obj3/paramB'] = -300
    obj_sys['obj3/paramB'] = -300

    print(obj_sys._reference_dict)

    pg = pgraph.ParameterGraph(sys)
    print(pg.object_path_built)

    print('----------')
    pg.print_parameters_eval()
    #print(pg.dict_parameters_eval())
    assert(pg.dict_parameters_eval() == {
        'o/p1': 100,
        'o/obj1/obj3/paramB': -300,
        'o/obj1/obj3/paramA': 300,
    })


def test_obj_coverage3():
    sys = pgraph.ParameterObject()
    sys['o/'] = obj_sys = pgraph.ParameterObject()
    obj_sys['obj1/'] = pgraph.ParameterObject()
    obj_sys['obj2/'] = pgraph.ParameterObject()
    obj_sys['obj1/obj3/'] = pgraph.ParameterObject()
    obj_sys['obj1/obj3/obj4/'] = pgraph.ParameterObject()
    #obj_sys.set_assign('obj3/', 'obj1/obj3/')
    obj_sys['obj3/'] = 'obj1/obj3/'

    obj_sys['p1'] = 100
    sys['o/obj3/obj4/test'] = 1

    obj1 = obj_sys['obj1/']
    obj1['imed/'] = pgraph.ParameterObject()
    obj1['obj3/paramA'] = 300
    obj2 = obj_sys['obj2/']

    #obj2['obj3/paramB'] = -300
    obj_sys['obj3/paramB'] = -300

    print(obj_sys._reference_dict)

    pg = pgraph.ParameterGraph(sys)
    print(pg.object_path_built)

    print('----------')
    pg.print_parameters_eval()
    #print(pg.dict_parameters_eval())
    assert(pg.dict_parameters_eval() == {
        'o/p1': 100,
        'o/obj1/obj3/paramB': -300,
        'o/obj1/obj3/paramA': 300,
        'o/obj1/obj3/obj4/test': 1,
    })
