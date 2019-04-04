#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from manage_crystal import file_parser

THIS_DIR = os.path.dirname(__file__)

def test_cif_parser():
    # test file with uncertainty
    crystal = file_parser.parse_from_filepath(os.path.join('input_cif', 'LEV.cif'), tm=None)
    assert crystal.length[0] == 13.1680
    assert crystal.length[1] == 13.1680
    assert crystal.length[2] == 22.5780
    assert crystal.angle_deg[0] == 90.0000
    assert crystal.angle_deg[1] == 90.0000
    assert crystal.angle_deg[2] == 120.0000

    # test file without uncertainty
    crystal = file_parser.parse_from_filepath(os.path.join('input_cif', 'core-cof.cif'), tm=None)
    assert crystal.length[0] == 23.4486008
    assert crystal.length[1] == 23.4486008
    assert crystal.length[2] == 3.9693000
    assert crystal.angle_deg[0] == 90.0000000
    assert crystal.angle_deg[1] == 90.0000
    assert crystal.angle_deg[2] == 120.0000
