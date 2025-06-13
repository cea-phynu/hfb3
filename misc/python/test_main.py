# ==============================================================================
#    HFB3
#    Copyright CEA, DAM F-91297 Arpajon, France
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.    If not, see <http://www.gnu.org/licenses/>.
# ==============================================================================

import hfb3

# ==============================================================================


def noel_test(func):
    def wrapper():
        print(f".. testing {func.__name__}")
        func()
        print(f"OK testing {func.__name__}")
    return wrapper

# ==============================================================================
# ==============================================================================
# ==============================================================================


def test_energy_contributions():
    filename = "examples/42Ca_deformed_2x9.msg.gz"

    dataTree = hfb3.DataTree(filename)
    state = hfb3.State(dataTree)
    interaction = hfb3.Interaction(dataTree, state)
    interaction.calcEnergies()

    err = 1e-6

    # print(f'{interaction("kinetic"         ).energy(hfb3.NEUTRON, hfb3.Field.DIRECT  )}')
    # print(f'{interaction("coulomb (Slater)").energy(hfb3.NEUTRON, hfb3.Field.DIRECT  )}')
    # print(f'{interaction("central", 0      ).energy(hfb3.NEUTRON, hfb3.Field.DIRECT  )}')
    # print(f'{interaction("central", 1      ).energy(hfb3.NEUTRON, hfb3.Field.DIRECT  )}')
    # print(f'{interaction("central", 0      ).energy(hfb3.NEUTRON, hfb3.Field.EXCHANGE)}')
    # print(f'{interaction("central", 1      ).energy(hfb3.NEUTRON, hfb3.Field.EXCHANGE)}')
    # print(f'{interaction("spin-orbit"      ).energy(hfb3.NEUTRON, hfb3.Field.DIRECT  )}')
    # print(f'{interaction("density"         ).energy(hfb3.NEUTRON, hfb3.Field.DIRECT  )}')
    # print(f'{interaction("2-body COM cor." ).energy(hfb3.NEUTRON, hfb3.Field.DIRECT  )}')
    # print(f'{interaction("2-body COM cor." ).energy(hfb3.NEUTRON, hfb3.Field.PAIRING )}')
    # print(f'{interaction("coulomb (Slater)").energy(hfb3.NEUTRON, hfb3.Field.EXCHANGE)}')
    # print(f'{interaction("central", 0      ).energy(hfb3.NEUTRON, hfb3.Field.PAIRING )}')
    # print(f'{interaction("central", 1      ).energy(hfb3.NEUTRON, hfb3.Field.PAIRING )}')
    # print(f'{interaction("rearrangement"   ).energy(hfb3.NEUTRON, hfb3.Field.DIRECT  )}')
    # print(f'{interaction("kinetic"         ).energy(hfb3.PROTON , hfb3.Field.DIRECT  )}')
    # print(f'{interaction("coulomb (Slater)").energy(hfb3.PROTON , hfb3.Field.DIRECT  )}')
    # print(f'{interaction("central", 0      ).energy(hfb3.PROTON , hfb3.Field.DIRECT  )}')
    # print(f'{interaction("central", 1      ).energy(hfb3.PROTON , hfb3.Field.DIRECT  )}')
    # print(f'{interaction("central", 0      ).energy(hfb3.PROTON , hfb3.Field.EXCHANGE)}')
    # print(f'{interaction("central", 1      ).energy(hfb3.PROTON , hfb3.Field.EXCHANGE)}')
    # print(f'{interaction("spin-orbit"      ).energy(hfb3.PROTON , hfb3.Field.DIRECT  )}')
    # print(f'{interaction("density"         ).energy(hfb3.PROTON , hfb3.Field.DIRECT  )}')
    # print(f'{interaction("2-body COM cor." ).energy(hfb3.PROTON , hfb3.Field.DIRECT  )}')
    # print(f'{interaction("2-body COM cor." ).energy(hfb3.PROTON , hfb3.Field.PAIRING )}')
    # print(f'{interaction("coulomb (Slater)").energy(hfb3.PROTON , hfb3.Field.EXCHANGE)}')
    # print(f'{interaction("central", 0      ).energy(hfb3.PROTON , hfb3.Field.PAIRING )}')
    # print(f'{interaction("central", 1      ).energy(hfb3.PROTON , hfb3.Field.PAIRING )}')
    # print(f'{interaction("rearrangement"   ).energy(hfb3.PROTON , hfb3.Field.DIRECT  )}')

    assert abs(interaction("kinetic"         ).energy(hfb3.NEUTRON, hfb3.Field.DIRECT  ) - 370.7063672673921     ) < err
    assert abs(interaction("coulomb (Slater)").energy(hfb3.NEUTRON, hfb3.Field.DIRECT  ) - 0.0                   ) < err
    assert abs(interaction("central", 0      ).energy(hfb3.NEUTRON, hfb3.Field.DIRECT  ) - -1007.1241532035914   ) < err
    assert abs(interaction("central", 1      ).energy(hfb3.NEUTRON, hfb3.Field.DIRECT  ) - -44.992952058861974   ) < err
    assert abs(interaction("central", 0      ).energy(hfb3.NEUTRON, hfb3.Field.EXCHANGE) - 458.1643288145652     ) < err
    assert abs(interaction("central", 1      ).energy(hfb3.NEUTRON, hfb3.Field.EXCHANGE) - -529.4984207720548    ) < err
    assert abs(interaction("spin-orbit"      ).energy(hfb3.NEUTRON, hfb3.Field.DIRECT  ) - -7.154443599388073    ) < err
    assert abs(interaction("density"         ).energy(hfb3.NEUTRON, hfb3.Field.DIRECT  ) - 557.0563434138456     ) < err
    assert abs(interaction("2-body COM cor." ).energy(hfb3.NEUTRON, hfb3.Field.DIRECT  ) - 4.629347682559987     ) < err
    assert abs(interaction("2-body COM cor." ).energy(hfb3.NEUTRON, hfb3.Field.PAIRING ) - 0.20174258670395293   ) < err
    assert abs(interaction("coulomb (Slater)").energy(hfb3.NEUTRON, hfb3.Field.EXCHANGE) - 0.0                   ) < err
    assert abs(interaction("central", 0      ).energy(hfb3.NEUTRON, hfb3.Field.PAIRING ) - 4.329007759557796     ) < err
    assert abs(interaction("central", 1      ).energy(hfb3.NEUTRON, hfb3.Field.PAIRING ) - -8.604895883359408    ) < err
    assert abs(interaction("rearrangement"   ).energy(hfb3.NEUTRON, hfb3.Field.DIRECT  ) - 96.94064120019276     ) < err
    assert abs(interaction("kinetic"         ).energy(hfb3.PROTON , hfb3.Field.DIRECT  ) - 310.3215733441709     ) < err
    assert abs(interaction("coulomb (Slater)").energy(hfb3.PROTON , hfb3.Field.DIRECT  ) - 79.17077057548208     ) < err
    assert abs(interaction("central", 0      ).energy(hfb3.PROTON , hfb3.Field.DIRECT  ) - -1014.983196657679    ) < err
    assert abs(interaction("central", 1      ).energy(hfb3.PROTON , hfb3.Field.DIRECT  ) - -21.069868044009247   ) < err
    assert abs(interaction("central", 0      ).energy(hfb3.PROTON , hfb3.Field.EXCHANGE) - 451.1226385016961     ) < err
    assert abs(interaction("central", 1      ).energy(hfb3.PROTON , hfb3.Field.EXCHANGE) - -516.2029836301575    ) < err
    assert abs(interaction("spin-orbit"      ).energy(hfb3.PROTON , hfb3.Field.DIRECT  ) - -1.6851980874629398   ) < err
    assert abs(interaction("density"         ).energy(hfb3.PROTON , hfb3.Field.DIRECT  ) - 557.0563454575162     ) < err
    assert abs(interaction("2-body COM cor." ).energy(hfb3.PROTON , hfb3.Field.DIRECT  ) - 3.6782755305554016    ) < err
    assert abs(interaction("2-body COM cor." ).energy(hfb3.PROTON , hfb3.Field.PAIRING ) - 1.9522566366369246e-13) < err
    assert abs(interaction("coulomb (Slater)").energy(hfb3.PROTON , hfb3.Field.EXCHANGE) - -7.449355269633595    ) < err
    assert abs(interaction("central", 0      ).energy(hfb3.PROTON , hfb3.Field.PAIRING ) - 3.799943157779036e-12 ) < err
    assert abs(interaction("central", 1      ).energy(hfb3.PROTON , hfb3.Field.PAIRING ) - -7.681252310532097e-12) < err
    assert abs(interaction("rearrangement"   ).energy(hfb3.PROTON , hfb3.Field.DIRECT  ) - 88.74405645476861     ) < err

# ==============================================================================
# ==============================================================================
# ==============================================================================


@noel_test
def test_fragment_properties():
    filename = "examples/42Ca_deformed_2x9.msg.gz"

    dataTree = hfb3.DataTree(filename)

    solverHFBBroyden = hfb3.SolverHFBBroyden(dataTree)
    discrete = hfb3.Discrete(solverHFBBroyden.state.basis,
                             hfb3.Mesh.regular(0, 0, -20, 10, 0, 20, 101, 1, 201))
    denst = discrete.getLocalXZ(solverHFBBroyden.state.rho(hfb3.NEUTRON) + solverHFBBroyden.state.rho(hfb3.PROTON), True)
    geomt = hfb3.Geometry(discrete.mesh, denst, hfb3.System(dataTree))
    izNeck = geomt.izNeck

    # if no neck is found, split at z = 0
    if izNeck < 0:
        izNeck = 100

    geomf = hfb3.Geometry(discrete.mesh, denst, hfb3.System(dataTree), izNeck)

    assert abs(geomf.intLeft                  - 21.7333865804) < 1e-2
    assert abs(geomf.intRight                 - 20.2592710507) < 1e-2
    assert abs(geomf.intLeft + geomf.intRight - 41.9926576312) < 2e-2

# ==============================================================================
# ==============================================================================
# ==============================================================================


@noel_test
def test_datatree():

    vq20 = 999
    myDict = {'constraints/q20t': vq20}

    # DataTree merging. The DataTree on the right of the '+' sign may overwrite values
    d = hfb3.DataTree("examples/42Ca_deformed_2x9.msg.gz") + hfb3.dictToDataTree(myDict)

    assert not d.contains('tutu')
    assert d.contains('system/nProt')
    assert d.getD('constraints/q10t') == 0
    assert d.getD('constraints/q20t') == vq20

# ==============================================================================
# ==============================================================================
# ==============================================================================


@noel_test
def test_dictToDataTree():
    import numpy as np
    from numpy.linalg import norm

    d = hfb3.dictToDataTree({
        'int'  : 1,
        'dbl'  : 2.0,
        'str'  : '3.00',
        'vec'  : np.ones((3,     ), dtype=np.float64, order='F'),
        'mat'  : np.ones((3, 3   ), dtype=np.float64, order='F'),
        'cube' : np.ones((3, 3, 3), dtype=np.float64, order='F'),
        'ivec' : np.ones((3,     ), dtype=np.int32  , order='F'),
        'imat' : np.ones((3, 3   ), dtype=np.int32  , order='F'),
        'icube': np.ones((3, 3, 3), dtype=np.int32  , order='F'),
        'uvec' : np.ones((3,     ), dtype=np.uint64 , order='F'),
        'umat' : np.ones((3, 3   ), dtype=np.uint64 , order='F'),
        'ucube': np.ones((3, 3, 3), dtype=np.uint64 , order='F'), })

    assert d.getI('int') == 1
    assert d.getD('dbl') == 2.0
    assert d.getS('str') == '3.00'
    assert norm(d.getV('vec'   ) - np.ones((3,     ), dtype=np.float64, order='F')) < 1e-14
    assert norm(d.getM('mat'   ) - np.ones((3, 3   ), dtype=np.float64, order='F')) < 1e-14
    assert norm(d.getC('cube'  ) - np.ones((3, 3, 3), dtype=np.float64, order='F')) < 1e-14
    assert norm(d.getIV('ivec' ) - np.ones((3,     ), dtype=np.int32  , order='F')) < 1e-14
    assert norm(d.getIM('imat' ) - np.ones((3, 3   ), dtype=np.int32  , order='F')) < 1e-14
    assert norm(d.getIC('icube') - np.ones((3, 3, 3), dtype=np.int32  , order='F')) < 1e-14
    assert norm(d.getUV('uvec' ) - np.ones((3,     ), dtype=np.uint64 , order='F')) < 1e-14
    assert norm(d.getUM('umat' ) - np.ones((3, 3   ), dtype=np.uint64 , order='F')) < 1e-14
    assert norm(d.getUC('ucube') - np.ones((3, 3, 3), dtype=np.uint64 , order='F')) < 1e-14

# ==============================================================================
# ==============================================================================
# ==============================================================================


@noel_test
def test_zippedDataTree():
    import numpy as np
    import gzip

    d = hfb3.dictToDataTree({
        'int'  : 1,
        'dbl'  : 2.0,
        'str'  : '3.00',
        'vec'  : np.ones((3,     ), dtype=np.float64, order='F'),
        'mat'  : np.ones((3, 3   ), dtype=np.float64, order='F'),
        'cube' : np.ones((3, 3, 3), dtype=np.float64, order='F'),
        'ivec' : np.ones((3,     ), dtype=np.int32  , order='F'),
        'imat' : np.ones((3, 3   ), dtype=np.int32  , order='F'),
        'icube': np.ones((3, 3, 3), dtype=np.int32  , order='F'), })

    serialized = hfb3.dataTreeToBytes(d)
    assert type(serialized) is bytes

    contentz = gzip.compress(serialized)
    assert type(contentz) is bytes

    d2 = hfb3.bytesToDataTree(serialized)
    assert d == d2

    content = gzip.decompress(contentz)
    d3 = hfb3.bytesToDataTree(content)
    assert d == d3

    d4 = hfb3.bytesToDataTree(contentz)
    assert d == d4

# ==============================================================================
# ==============================================================================
# ==============================================================================


@noel_test
def test_dataTreeFromZippedContent():
    filename = "examples/42Ca_deformed_2x9.msg.gz"
    with open(filename, 'rb') as f:
        content = f.read()
        dataTree = hfb3.bytesToDataTree(content)
        assert dataTree.getI('basis/nOscil') == 9

# ==============================================================================
# ==============================================================================
# ==============================================================================


@noel_test
def test_dataTreeFromUnzippedContent():
    import gzip
    filename = "examples/42Ca_deformed_2x9.msg.gz"
    with gzip.open(filename, 'rb') as f:
        content = f.read()
        dataTree = hfb3.bytesToDataTree(content)
        assert dataTree.getI('basis/nOscil') == 9

# ==============================================================================
# ==============================================================================
# ==============================================================================


@noel_test
def test_state():
    filename = "examples/42Ca_deformed_2x9.msg.gz"

    # test file -> constraints
    dataTree = hfb3.DataTree(filename)
    state = hfb3.State(dataTree)
    assert state.constraints['q10t'].val ==    0.0
    assert state.constraints['q20t'].val == 10.0

    # test changing constraint values
    state.constraints['q10t'].val =    40.0
    state.constraints['q20t'].val = -23.0
    assert state.constraints['q10t'].val ==    40.0
    assert state.constraints['q20t'].val == -23.0

    # test changing constraint values from a Python dict
    state.constraints = hfb3.map_constraints({'q10t': hfb3.Constraint('q10t', 42.0),
                                                        'q20t': hfb3.Constraint('q20t', 43.0)})
    assert state.constraints['q10t'].val ==    42.0
    assert state.constraints['q20t'].val ==    43.0

    # test removing a constraint
    dataTree2 = hfb3.dictToDataTree({'constraints/q20t': None})
    dataTree.merge(dataTree2)
    state = hfb3.State(dataTree)
    assert state.constraints.keys() == ['q10t']


# ==============================================================================
# ==============================================================================
# ==============================================================================


@noel_test
def test_qnumbers():
    import numpy as np

    qn = hfb3.Qnumbers(2)
    qn.setNames(('a,', 'b'))
    qn.append(np.array([0, 0], dtype=np.dtype('i4')))
    qn.append(np.array([0, 1], dtype=np.dtype('i4')))
    qn.append(np.array([1, 0], dtype=np.dtype('i4')))
    qn.append(np.array([1, 1], dtype=np.dtype('i4')))

    qn.calcBlocks(np.array([1,], dtype=np.dtype('u8')))
    assert qn.blocks.size() == 2
    assert qn.blocks[0].nb == 2
    assert np.all(qn.blocks[0](0) == np.array([0, 0]))
    assert np.all(qn.blocks[0](1) == np.array([1, 0]))
    assert np.all(qn.blocks[1](0) == np.array([0, 1]))
    assert np.all(qn.blocks[1](1) == np.array([1, 1]))

# ==============================================================================
# ==============================================================================
# ==============================================================================


@noel_test
def test_multi():
    import numpy as np
    from numpy.linalg import norm

    multiInteger = hfb3.MultiI()
    multiInteger[(2, -6, 0, 0, 87)] = -1234
    assert multiInteger[(2, -6, 0, 0, 87)] == multiInteger(2, -6, 0, 0, 87)

    multiDouble = hfb3.MultiD()
    multiDouble[(3, 8, 1)] = 9.87
    multiDouble[(5, 5   )] = 6.2
    multiDouble2 = hfb3.MultiD()
    multiDouble2[(5, 5   )] = 6.2
    multiDouble2[(3, 8, 1)] = 9.87
    assert multiDouble.info() == multiDouble2.info()
    assert multiDouble == multiDouble2

    # types V, IV, UV

    multiV = hfb3.MultiV()
    multiV[(2, -6, 87)] = np.zeros((4,), dtype=np.float64, order="F")
    assert norm(multiV[(2, -6, 87)] - multiV(2, -6, 87)) < 1e-14

    multiIV = hfb3.MultiIV()
    multiIV[(2, -6, 87)] = np.zeros((4,), dtype=np.int32, order="F")
    multiIV[(2, -2)] = np.zeros((0,), dtype=np.int32, order="F")
    assert multiIV.size() == 2

    multiUV = hfb3.MultiUV()
    multiUV[(2, -6, 87)] = np.zeros((4,), dtype=np.uint64, order="F")
    multiUV[(2, -2)] = np.zeros((0,), dtype=np.uint64, order="F")
    assert multiUV.size() == 2

    # types M, IM, UM

    multiM = hfb3.MultiM()
    multiM[(2, -6, 87)] = np.zeros((4, 3), dtype=np.float64, order="F")
    assert norm(multiM[(2, -6, 87)] - multiM(2, -6, 87)) < 1e-14

    multiIM = hfb3.MultiIM()
    multiIM[(2, -6, 87)] = np.zeros((4, 2), dtype=np.int32, order="F")
    multiIM[(2, -2)] = np.zeros((0, 0), dtype=np.int32, order="F")
    assert multiIM.size() == 2

    multiUM = hfb3.MultiUM()
    multiUM[(2, -6, 87)] = np.zeros((4, 2), dtype=np.uint64, order="F")
    multiUM[(2, -2)] = np.zeros((0, 0), dtype=np.uint64, order="F")
    assert multiUM.size() == 2

    # types C, IC, UC

    multiC = hfb3.MultiC()
    multiC[(2, -6, 87)] = np.zeros((4, 3, 2), dtype=np.float64, order="F")
    assert norm(multiC[(2, -6, 87)] - multiC(2, -6, 87)) < 1e-14

    multiIC = hfb3.MultiIC()
    multiIC[(2, -6, 87)] = np.zeros((4, 2, 3), dtype=np.int32, order="F")
    multiIC[(2, -2)] = np.zeros((0, 0, 0), dtype=np.int32, order="F")
    assert multiIC.size() == 2

    multiUC = hfb3.MultiUC()
    multiUC[(2, -6, 87)] = np.zeros((4, 2, 3), dtype=np.uint64, order="F")
    multiUC[(2, -2)] = np.zeros((0, 0, 0), dtype=np.uint64, order="F")
    assert multiUC.size() == 2

# ==============================================================================
# ==============================================================================
# ==============================================================================


@noel_test
def test_multiStates():

    mi0 = hfb3.MultiI()
    mi0[(5,)] = 8

    mi1 = hfb3.MultiI()
    mi1[(5,)] = 8

    assert mi0 == mi1

    ms0 = hfb3.MultiStates()
    ms0[(hfb3.NEUTRON,)] = hfb3.States()
    ms0[(hfb3.PROTON, )] = hfb3.States()

    ms1 = hfb3.MultiStates()
    ms1[(hfb3.NEUTRON,)] = hfb3.States()
    ms1[(hfb3.PROTON, )] = hfb3.States()

    assert ms0 == ms1

# ==============================================================================
# ==============================================================================
# ==============================================================================


if __name__ == '__main__':

    hfb3.cvar.useColors = False
    hfb3.cvar.msgToOut = [hfb3.MSG_ERROR]

    test_energy_contributions()
    test_fragment_properties()
    test_datatree()
    test_dictToDataTree()
    test_zippedDataTree()
    test_dataTreeFromZippedContent()
    test_dataTreeFromUnzippedContent()
    test_state()
    test_qnumbers()
    test_multi()
    test_multiStates()
