import eval_util

def compute_module_count(inf, is_mmod):
    if is_mmod:
        mods = eval_util.read_mmods(inf)
    else:
        mods = eval_util.read_rmods(inf)

    return len(mods)