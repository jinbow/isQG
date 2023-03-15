def invertit(filein: str, fileout: str) -> None:
    """
    This function reads in a MATLAB file containing ISQG data, solves the
    SQG and PSI equations using the isqg_data class, and saves the results to
    a new MATLAB file.

    Args:
        filein (str): Path to the input MATLAB file.
        fileout (str): Path to the output MATLAB file.

    Returns:
        None
    """
    from isqg import isqg_data
    from scipy.io import loadmat, savemat

    # Initialize isqg_data object
    d = isqg_data()

    # Load input MATLAB file
    f = loadmat(filein, squeeze_me=True)

    # Set attributes of isqg_data object to values from input file
    for key in f.keys():
        if key[:2] != '__':
            print("read", key)
            setattr(d, key, f[key])

    # Solve SQG and PSI equations using isqg_data object
    d.solve_sqg()
    d.solve_psii()

    # Save output variables to a new MATLAB file
    outvar = ['us', 'vs', 'ui', 'vi', 'ut', 'vt', 'psis', 'psii', 'psit',
              'ssh', 'ssd', 'rhos', 'rhoi', 'rho0', 'f0', 'z', 'zf', 'zc',
              'lon', 'lat', 'vorticity']

    dout = {}
    for ov in outvar:
        dout[ov] = getattr(d, ov)

    savemat(fileout, dout)


if __name__ == '__main__':
    import sys

    if len(sys.argv)<3:
        print("Usage: python invert.py filein.mat fileout.mat")
        sys.exit(0)
    else:
        filein,fileout = sys.argv[1:3]

    invertit(filein,fileout)
