

def binary_search_range(arr, target):
    # search for starting and final index of elements in neighbor list,
    # which has elements of [ idx1, idx2 ], where idx1 is the origin atom,
    # and idx2 is neighbor. The neighbor list has to be ordered by idx1 on entry.
    low = 0
    high = len(arr) - 1
    start_index = -1
    end_index = -1

    # Find the start index
    while low <= high:
        mid = (low + high) // 2
        if arr[mid][0] == target:
            start_index = mid
            high = mid - 1
        elif arr[mid][0] < target:
            low = mid + 1
        else:
            high = mid - 1

    low, high = 0, len(arr) - 1

    # Find the end index
    while low <= high:
        mid = (low + high) // 2
        if arr[mid][0] == target:
            end_index = mid
            low = mid + 1
        elif arr[mid][0] < target:
            low = mid + 1
        else:
            high = mid - 1

    return start_index, end_index

def extract_elements(neiglist, veclist, target):
    import numpy as np
    ## find the starting and final index of arr where first element == target
    start_index, end_index = binary_search_range(neiglist, target)

    res_idx = np.empty([0], dtype=int)
    res_vec = np.ndarray([0,0], dtype=np.float64)
    ## put the target idx into first place
    # res_idx.append( target )
    res_idx = np.append( res_idx, target )
    res_vec = np.append( res_vec, np.array([0.0,0.0,0.0]))


    ## extend the result with array values
    # res_idx.extend( neiglist[start_index:end_index+1, 1] )
    res_idx = np.append( res_idx, neiglist[start_index:end_index+1, 1])
    for v in veclist[start_index:end_index+1]:
        res_vec = np.vstack( [res_vec, v] )
    return res_idx, res_vec

