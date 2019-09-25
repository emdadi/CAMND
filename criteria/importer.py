import libsbml
import sys
import numpy
import my_constants
import sys
import math

if my_constants.win:
    import scipy.io


def __enoughInStr(z):
    return str(z if (type(z) is str or type(z) is unicode or type(z) is numpy.unicode_) \
                   else z[0] if (type(z[0]) is str or type(z[0]) is unicode or type(z[0]) is numpy.unicode_) \
        else z[0][0])


def __enoughInInt(z):
    return int(z if type(z) is int or type(z) is numpy.uint8 \
                   else z[0] if type(z[0]) is int or type(z[0]) is numpy.uint8 \
        else z[0][0])


# no more than 9 zeros between two fractional digits is acceptable! i.e. fract_len(1e-11) is zero but fract_len(1e-10) is 10!
def fract_len(f):
    l = 0
    while True:
        if l > 100:
            raise Exception('fract_len failed :(')

        z = f * 10 ** l
        if abs(math.modf(z)[0]) < 1e-10 or abs(math.modf(z)[0]) > 1 - 1e-10:
            return l

        l += 1


def normalize_list_of_floats_to_ints(fl):
    non_ints = []
    min_non_int = sys.maxsize
    most_fractional_non_int = -1
    for f in fl:
        if f != int(f):
            non_ints.append(f)
            if abs(f) < abs(min_non_int):
                min_non_int = f
            if fract_len(f) > fract_len(most_fractional_non_int):
                most_fractional_non_int = f

    ret = fl
    if non_ints:
        factor = 10 ** fract_len(most_fractional_non_int)
        for i in range(len(fl)):
            ret[i] = math.trunc(round(fl[i] * factor))
    else:
        ret = [math.trunc(fl[i]) for i in range(len(fl))]

    return ret


def convertToAllIntsEachReactIndependent(s, mc, rc):
    for ri in range(rc):
        r = [s[i][ri] for i in range(mc)]
        r = normalize_list_of_floats_to_ints(r)
        for i in range(mc):
            s[i][ri] = r[i]

    return s


def convertToAllIntsOld(s, mc, rc):
    nonints = []
    minni = (sys.maxint, -1, -1)
    for r in range(mc):
        for c in range(rc):
            if s[r][c] != int(s[r][c]):
                nonints.append((s[r][c], r, c))
                if abs(s[r][c]) < abs(minni[0]) and round(s[r][c]) != 0:
                    minni = (s[r][c], r, c)

    if minni[1] != -1:
        factor = 100.0 / round(abs(minni[0]))
        for i in range(len(s)):
            for j in range(len(s[0])):
                s[i][j] = int(round(s[i][j] * factor))

    return s


def matlabStiochiometricMatrix(srcFile, sIdx=0, revIdx=1, rxnsIdx=2, metsIdx=3, doPrint=False):
    # return S, met, rxns, revs
    m = scipy.io.loadmat(srcFile)
    key = filter(lambda x: not (x.startswith('__') or x.endswith('__')), m.keys())[0]
    mats = m[key][0][0]
    rxns = [__enoughInStr(z) for z in mats[rxnsIdx]]
    rc = len(rxns)
    mets = [__enoughInStr(z) for z in mats[metsIdx]]
    mc = len(mets)
    S = mats[sIdx]
    rev = [__enoughInInt(z) for z in mats[revIdx]]
    s = [[0] * rc for i in range(mc)]
    for r in range(mc):
        for c in range(rc):
            s[r][c] = S[r, c]

    if doPrint:
        print s

    s = convertToAllIntsEachReactIndependent(s, mc, rc)

    return s, mets, rxns, rev


def inlineStiochiometricMatrix(doPrint):
    # return S, metaboliteCount, reactionCount, externalMetaboliteCount
    #    cols = 9
    #    scoeffsmat = [
    #        1, -1,  0, -1,  0,  0,  0,  0,  0,
    #        0,  1, -1,  0,  1,  0,  0,  0,  0,
    #        0,  0,  0,  1, -1,  0,  0,  0,  0,
    #        0,  0,  1,  0,  0, -1,  0, -1,  0,
    #        0,  0,  0,  0,  0,  1, -1,  0,  1,
    #        0,  0,  0,  0,  0,  0,  0,  1, -1,
    #        ]

    cols = 11
    scoeffsmat = [1, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 1,

                  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, ]

    #    cols = 11
    #    scoeffsmat = [
    #        3,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,
    #        1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    #        0,  0, -1,  0, -1,  0,  0, -1,  0,  1,  0,
    #        0,  1,  0, -1,  0,  0,  0,  0,  0,  0,  0,
    #        0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  0,
    #        0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,
    #        0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,
    #        0,  0,  0,  2,  0,  0,  1,  0,  0,  0, -1,
    #        0,  0,  0,  0,  0,  0,  0,  0, -3,  0,  1,
    #        0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,
    #    ]

    reactionCount = cols
    metaboliteCount = len(scoeffsmat) / cols
    externalMetaboliteCount = 2

    S = lp.pylinpro.matlib.matzero(metaboliteCount, reactionCount)

    lp.pylinpro.matlib.matsetblock(S, 0, 0, metaboliteCount - 1, reactionCount - 1, scoeffsmat)

    if doPrint:
        lp.pylinpro.matlib.matprint(S, "%4d", [chr(i) for i in range(ord('A'), ord('Z'))], range(1, reactionCount - 1) + ['I', 'J'], "%4s")

    return S, [chr(ord('A') + x) for x in xrange(metaboliteCount)], [str(x + 1) for x in xrange(reactionCount)], [0] * reactionCount

    S = [[0, 0, 0, 0, -100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -100, 0, 100, 100, -100, 0, -100, 0, 0, 0, -100, 0, 0, -100, -100, 0, -100, 100, 0, 0, 0, 0, 0, -100, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, -100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 100, -100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, -100, 0, 0, 0, 0, 0],
         [0, 0, -100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, -100, 0, 100, 0, 0, 0, 0, 0, 100, 0, 0, 100, 0, 0, 0, -100, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -100, -100, 100, 0, 0, 0, 0, 0, 0, -100, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 200, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -100, -100, -100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 200, 200, -300, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3],
         [-100, -100, 0, -100, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 100, 100, 0, 0, 100, 0, 0, 0, 0, 0, 0, 100, -100, -100, -100, -100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 100, -100, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [-100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -100, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -100, -100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]

    mets = ['2OG_p', '2PG_c', '3PG_c', 'AcSer_p', 'CoA_p', 'F6P_p', 'Gaba_c', 'Glc_p', 'H_ext', 'Ile_p', 'Mal_c', 'Ser_c', 'Succ_p', 'THDPA_p']

    rxns = ['R1007', 'R1008', 'R1009', 'R1010', 'R1012', 'R1014', 'R1015', 'R1016', 'R1017', 'R1021', 'R742', 'R743', 'R744', 'R747', 'R748', 'R749', 'R750', 'R751', 'R752', 'R755', 'R756', 'R757', 'R758', 'R759', 'R760', 'R761', 'R763', 'R764', 'R769', 'R774', 'R777', 'R780', 'R781', 'R782', 'R783', 'R794', 'R795', 'R803', 'R804', 'R805', 'R814', 'R816', 'R821', 'R822', 'R823', 'R825', 'R828', 'R830', 'R833', 'R837', 'R848', 'R849', 'R850', 'R852', 'R853', 'R854', 'R855', 'R857', 'R860',
            'R862', 'R865', 'R866', 'R867', 'R869', 'R870', 'R872', 'R873', 'R881', 'R884', 'R885', 'R886', 'R888', 'R894', 'R895', 'R896', 'R897', 'R907', 'R910', 'R917', 'R920', 'R924', 'R925', 'R928', 'R940', 'R942', 'R943', 'R944', 'R945', 'R946', 'R952', 'R956', 'R957', 'R960', 'R962', 'R965', 'R966', 'R967', 'R968', 'R973', 'R975', 'R998']

    revs = [1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1]

    S = [[100, -100, -100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-100, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
    mets = ['HOMO-Cys_c', 'SAH_c']
    rxns = ['R1003', 'R1005', 'R1006', 'R1015', 'R1016', 'R1017', 'R1020', 'R778', 'R789', 'R790', 'R792', 'R793', 'R843', 'R844', 'R846', 'R847', 'R848', 'R849', 'R851', 'R852', 'R856', 'R857', 'R879', 'R880', 'R892', 'R899', 'R904', 'R908', 'R909', 'R931', 'R944', 'R945', 'R946', 'R947', 'R949', 'R953', 'R956', 'R964', 'R967', 'R972', 'R976', 'R998']
    revs = [0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1]

    return S, mets, rxns, revs


def csvStiochiometricMatrix(srcFile, doPrint):
    # return S, metaboliteCount, reactionCount
    f = open(srcFile, 'r')

    s = []
    for l in f:
        if l.strip() == '':
            break
        s.append(l.split(';'))

    rev = []
    reactions = []
    metabolites = []
    for l in f:
        if l.lower().startswith('REV'):
            rev = l.split(';')
        elif l.lower().startswith('RXNS'):
            reactions = l.split(';')
        elif l.lower().startswith('METS'):
            metabolites = l.split(';')

    mc = len(metabolites)
    rc = len(reactions)

    if doPrint:
        print s

    nonints = []
    minni = (sys.maxint, -1, -1)
    for r in range(mc):
        for c in range(rc):
            if s[r][c] != int(s[r][c]):
                nonints.append((s[r][c], r, c))
                if abs(s[r][c]) < abs(minni[0]) and round(s[r][c]) != 0:
                    minni = (s[r][c], r, c)

    if minni[1] != -1:
        factor = 100.0 / round(abs(minni[0]))
        for i in range(len(s)):
            for j in range(len(s[0])):
                s[i][j] = int(round(s[i][j] * factor))

    return s, metabolites, reactions, rev


def sbmlStoichiometricMatrix(srcFile, do_print=False, read_species_compart=False, remove_biomass=True, normalize_stoich=True):
    # return S, met, rxns, revs

    reader = libsbml.SBMLReader()
    doc = reader.readSBML(srcFile)
    m = doc.getModel()

    if do_print:
        for i in xrange(doc.getNumErrors()):
            if doc.getError(i).getErrorId() == libsbml.XMLFileUnreadable:
                print 'SBML: XMLFileUnreadable'
            elif doc.getError(i).getErrorId() == libsbml.XMLFileOperationError:
                print 'SBML: XMLFileOperationError'
            else:
                print 'SBML: ' + str(doc.getError(i).getErrorId())

        print 'SBML: numSpecies: %d, numReactions: %d' % (m.getNumSpecies(), m.getNumReactions())

    ms = m.getListOfSpecies()
    rs = m.getListOfReactions()

    mets = sorted([s.getId() for s in ms])
    metNames = [m.getSpecies(met).getName() for met in mets]
    if read_species_compart:
        comparts = [m.getSpecies(met).getCompartment() for met in mets]
    rxns = sorted([r.getId() for r in rs])
    rxnNames = [m.getReaction(rxn).getName() for rxn in rxns]
    revs = [1 if m.getReaction(rid).getReversible() else 0 for rid in rxns]

    metsMap = dict((m, i) for (i, m) in enumerate(mets))  # {m: i for (i,m) in enumerate(mets)}
    rxnsMap = dict((r, i) for (i, r) in enumerate(rxns))  # {r: i for (i,r) in enumerate(rxns)}

    S = [[0 for j in xrange(len(rxns))] for i in xrange(len(mets))]

    for r in rs:
        ridx = rxnsMap[r.getId()]

        reacs = r.getListOfReactants()
        for re in reacs:
            midx = metsMap[re.getSpecies()]
            S[midx][ridx] = -re.getStoichiometry()

        prods = r.getListOfProducts()
        for pr in prods:
            midx = metsMap[pr.getSpecies()]
            S[midx][ridx] = pr.getStoichiometry()

    biomass = []

    sbml_f = open(srcFile, 'r')
    sbml_text = sbml_f.read()
    obj_idx = sbml_text.find('fbc:fluxObjective')
    if obj_idx != -1:
        b_name_stidx = sbml_text.index('fbc:reaction="', obj_idx) + len('fbc:reaction="')
        b_name_edidx = sbml_text.index('"', b_name_stidx)
        b_name = sbml_text[b_name_stidx: b_name_edidx]

        biomass = [b_name, [S[i][rxns.index(b_name)] for i in range(len(mets))]]

        obj_idx = sbml_text.find('fbc:fluxObjective', b_name_edidx)
        if obj_idx != -1:
            raise Exception('more than one biomass for %s' % srcFile)

        if remove_biomass:
            all_but_biomass = list(rxns)
            all_but_biomass.remove(b_name)
            S, mets, rxns, revs = reduce_network_to_reacts(all_but_biomass, S, mets, rxns, revs)

    sbml_f.close()

    if normalize_stoich:
        S = convertToAllIntsEachReactIndependent(S, len(mets), len(rxns))

    # n1 = -1
    # n2 = -1
    # mm = -1
    # for zzai, zza in enumerate(S):
    #     for zzbi, zzb in enumerate(zza):
    #         if abs(zzb) > mm:
    #             mm = abs(zzb)
    #             n1 = zzai
    #             n2 = zzbi
    # print mm
    # print [S[i][n2] for i in range(len(mets))]

    if read_species_compart:
        return S, mets, rxns, revs, metNames, rxnNames, biomass, comparts
    else:
        return S, mets, rxns, revs, metNames, rxnNames, biomass


def reduce_network_to_reacts(mod_reacts, S, mets, reacts, revs):
    react_filter = [True if reacts[i] in mod_reacts else False for i in range(len(reacts))]
    Sp = [[S[i][j] for j in range(len(reacts)) if react_filter[j]] for i in range(len(mets))]
    metsP = list(mets)
    reactsP = [reacts[i] for i in range(len(reacts)) if react_filter[i]]
    revsP = [revs[i] for i in range(len(revs)) if react_filter[i]]

    return Sp, metsP, reactsP, revsP


def reduce_network_to_metabs(mod_metabs, S, mets, reacts, revs):
    metab_filter = [True if mets[i] in mod_metabs else False for i in range(len(mets))]
    Sp = [[S[i][j] for j in range(len(reacts))] for i in range(len(mets)) if metab_filter[i]]
    metsP = [mets[i] for i in range(len(mets)) if metab_filter[i]]
    reactsP = list(reacts)
    revsP = list(revs)

    return Sp, metsP, reactsP, revsP
