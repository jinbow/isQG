#!/usr/bin/python

def invertit(filein, fileout):
    from isqg import isqg_data
    from scipy.io import loadmat, savemat
    
    
    d = isqg_data()
    
    f = loadmat(filein,squeeze_me=True)
    
    for key in f.keys():
        if key[:2]!='__':
            print "read ",key
            setattr(d,key,f[key])
    
    d.solve_sqg()
    d.solve_psii()
    
    
    outvar = ['us','vs','ui','vi','ut','vt','psis','psii','psit','ssh','ssd',
              'rhos','rhoi','rho0','f0','z','zf','zc',
              'lon','lat','vorticity']
    
    dout = {}
    for ov in outvar:
        dout[ov] = getattr(d,ov)
    
    savemat(fileout, dout)

if __name__ == '__main__':
    import sys
    
    if len(sys.argv)<3:
        print "Usage: python invert.py filein.mat fileout.mat"
        sys.exit(0)
    else:
        filein,fileout = sys.argv[1:3]
    
    invertit(filein,fileout)