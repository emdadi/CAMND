import os
import time

import eval_util
from util import my_constants
from util.my_util import write_line
from util import my_util


def compute_and_print_size_distribution(inf, is_mmod, species, method):
    out_dir = '%s/sizes' % my_constants.evalResultPath
    my_util.mkdir_p(out_dir)

    if is_mmod:
        mods = eval_util.read_mmods(inf)
    else:
        mods = eval_util.read_rmods(inf)

    sizes = []
    for mod_name, mod_metabs in mods.iteritems():
        sizes.append(len(mod_metabs))

    result_file = '%s/%s_%s.png' % (out_dir, species, method)
    alt_result_file = '%s/sc_%s_%s.png' % (out_dir, species, method)

    f = open('size_distrib.r', 'w')
    write_line(f, 'png("%s")' % result_file)
    write_line(f, "hist(c(%s), main='%d modules (%s)', xlab='', ylab='')" % (','.join([str(s) for s in sizes]), len(mods), species + ' ' + method))
    write_line(f, 'dev.off()')
    write_line(f, 'png("%s")' % alt_result_file)
    write_line(f, "stripchart(c(%s), main='%d modules (%s)', method='stack', offset=0.5, pch=1)" % (','.join([str(s) for s in sizes]), len(mods), species + ' ' + method))
    write_line(f, 'dev.off()')
    f.close()

    res = os.system('%s size_distrib.r' % my_constants.rScriptPath)

    while not os.path.isfile(result_file):
        time.sleep(5)

    return '=HYPERLINK("%s")' % result_file


if __name__ == '__main__':
    # species = 'toy_model'
    species = 'ecoli_core'

    out_dir = r'%s/%s/newman' % (my_constants.resultPath, species)

    print compute_and_print_size_distribution('%s/final_modules.txt' % out_dir, True)