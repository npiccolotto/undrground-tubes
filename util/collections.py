def merge_alternating(xs, ys=None):
    if ys is None and len(xs) == 2:
        xs, ys = xs
    # https://stackoverflow.com/a/3678938/490524
    result = [None] * (len(xs) + len(ys))
    result[::2] = xs
    result[1::2] = ys
    return result
