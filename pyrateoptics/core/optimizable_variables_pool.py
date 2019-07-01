#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

from .optimizable_variable import OptimizableVariable
from .functionobjects_pool import FunctionObjectsPool
from .log import BaseLogger


class OptimizableVariablesPool(BaseLogger):

    def __init__(self, variables_pool_from_serialization_iter, name=""):
        super(OptimizableVariablesPool, self).__init__(
                name=name,
                kind="optimizablevariablespool")
        self.variables_pool = variables_pool_from_serialization_iter
        self.variables_pool = self.extendPool(self.variables_pool)
        # guarantee that pool contains ALL variables which are referenced
        # by pickups

    def extendPool(self, pool):
        """
        Collects all pickup variables in pool and checks whether their
        arguments are also in pool. If this is not the case than it adds them.

        @param pool: dictionary of unique_id strings with optimizable variables

        """

        def are_all_args_in_pool(pool):
            """
            Returns the truth value of whether all args of pickups are in the
            pool. Also returns the pickup variables in pool.

            @param pool: dictionary of unique strings with optimizable vars
            """
            res = False
            all_pickup_vars = list([ov for (k, ov) in pool.items()
                                    if ov.var_type() == "pickup"])
            res = all([all(
                    [a.unique_id in pool
                        for a in ov._state.parameters["args"]]
                    ) for ov in all_pickup_vars])

            return (res, all_pickup_vars)

        result_pool = {}
        # update result pool to pool items
        for (key, ov) in pool.items():
            result_pool[key] = ov

        # check if all of the args are in pool
        (all_args_in_pool, pickup_vars) = are_all_args_in_pool(result_pool)
        # while this is not the case, add the args to the pool
        while not all_args_in_pool:
            for ov in pickup_vars:
                for arg_variable in ov._state.parameters["args"]:
                    if arg_variable.unique_id not in result_pool:
                        result_pool[arg_variable.unique_id] = arg_variable
            (all_args_in_pool, pickup_vars) = are_all_args_in_pool(result_pool)
        return result_pool

    def generateFunctionObjectsPool(self, name=""):
        functionobjects_dictionary = {}
        for (key, ov) in self.variables_pool.items():
            (tfo, ftrafo, finvtrafo) = ov._transform_functionobject_functionname_triple
            if tfo.unique_id not in functionobjects_dictionary:
                functionobjects_dictionary[tfo.unique_id] = tfo
            if ov.var_type() == "pickup":
                (pfo, functionname) = ov._state.parameters["functionobject"]
                if pfo.unique_id not in functionobjects_dictionary:
                    functionobjects_dictionary[pfo.unique_id] = pfo
        return FunctionObjectsPool(functionobjects_dictionary, name=name)

    def toDictionary(self):
        return dict([(uid, variable.toDictionary())
                    for (uid, variable) in self.variables_pool.items()])

    @staticmethod
    def fromDictionary(dictionary, functionpool_dictionary,
                       source_checked,
                       variables_checked, name=""):
        """
        Reconstruct variable pool from dictionary. Pre initialize pickups
        and successively construct them such that at the end of the day
        a valid pool of variables is created.

        @param dictionary - pool dictionary which contains the structural
        information to reconstruct the variables pool.

        @param source_checked - True or False if source code was inspected by
                                user and considered safe.
        @param variables_checked - True or False if variables were inspected
                                    by user and considered safe.

        """

        def are_all_args_opt_vars(pool):
            res = False
            all_pickup_vars = list([ov for (k, ov) in pool.items()
                                    if ov.var_type() == "pickup"])
            res = all([all(
                    [isinstance(a, OptimizableVariable) for a in ov._state.parameters["args"]]
                    ) for ov in all_pickup_vars])

            return (res, all_pickup_vars)

        functionobjects_pool = FunctionObjectsPool.fromDictionary(
            functionpool_dictionary, source_checked, variables_checked)

        variables_pool = {}
        for (uid, var_dict) in dictionary.items():
            variables_pool[uid] = OptimizableVariable.fromDictionary(
                    var_dict, functionobjects_pool)

        (all_args_opt_vars, pickupvars) = are_all_args_opt_vars(variables_pool)
        while not all_args_opt_vars:
            for ov in pickupvars:
                new_args = []
                for arg in ov._state.parameters["args"]:
                    if isinstance(arg, str):
                        new_args.append(variables_pool[arg])
                    else:
                        new_args.append(arg)
                ov._state.parameters["args"] = new_args
                (all_args_opt_vars, pickupvars) = are_all_args_opt_vars(variables_pool)

        for (k, ov) in variables_pool.items():
            if ov.var_type() == "pickup":
                ov._state.isvalid = True

        return OptimizableVariablesPool(variables_pool, name=name)
