# aux
def read_xyz(fname):
    from numpy import genfromtxt, ndarray, float64, int32
    import linecache
    typ = genfromtxt( fname, skip_header=2, usecols=[0], dtype=int32 )
    coords = genfromtxt( fname, skip_header = 2, usecols = [1,2,3], dtype=float64 )
    nat = len( typ )
    # cc=ndarray( (nat,3), dtype=float64, order="C" )
    # cc=coords
    opts = linecache.getline( fname, 2 )
    options = opts.split()
    return nat, typ, coords, options


def count_nn( neiglist ):
    import numpy as np
    # get count of neighbours for each atom
    u_val, cnt = np.unique( neiglist[:,0], return_counts=True )
    # perform cumulative sum
    cnt_sum = np.cumsum( cnt )
    return cnt_sum

def extract_elements( neiglist, veclist, nn_sum, target ):
    '''
    nn_sum is the cumulative sum of neighbors
    '''
    import numpy as np
    start_index = nn_sum[ target - 1 ]
    end_index = nn_sum[ target ]
    if target == 0:
        start_index = 0

    # total size +1 because first element is not counted
    n = end_index - start_index + 1

    # declate arrays
    res_idx = np.empty([n], dtype=np.int32)
    res_vec = np.ndarray([n,3], dtype=np.float64)

    # fill first vector as (0.0, 0.0, 0.0)
    res_vec[0] = np.array([0.0, 0.0, 0.0])
    # copy other vectors
    for i,v in enumerate( veclist[start_index:end_index] ):
        res_vec[i+1] = v

    # first index is target
    res_idx[0] = target
    # copy other indices
    res_idx[1:] = neiglist[start_index:end_index, 1]

    return res_idx, res_vec

