# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 15:42:24 2019

@author: luiggi
"""
#-----------------------------------------------------------
# PARA DEFINIR EL PATH ABSOLUTO DE LOS MÓDULOS DE PYNOXTLI
#
import os, sys
if not("pynoxtli/base" in sys.path[0][-13:]):
    sys.path.insert(0, os.path.abspath('..'))
#-----------------------------------------------------------

import la.csr_fvm_ijk_1D as CSR_FVM_ijk_1D
import la.csr_fvm_ijk_2D as CSR_FVM_ijk_2D
