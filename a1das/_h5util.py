"""
A few utilities used for Febus and reducted file to get/set chunk cache size
"""
def _get_h5_chunk_cache_size(fid):
    """
    fid = H5py high level file descriptor

    return the chunk cache size
    """
    cacheSettings = list(fid.id.get_access_plist().get_cache())

    return(cacheSettings[2])
