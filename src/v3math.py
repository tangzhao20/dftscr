# Here are some basic operations for vectors

def v3dpp(x1, x2) :
    out=0.0
    for ix in range(3) :
        out+=(x1[ix]-x2[ix])*(x1[ix]-x2[ix])
    return out**0.5

def v3matchpp(x1, x2, tol=1e-6) :
    for ix in range(3) :
        if abs(x1[ix]-x2[ix])>tol :
            return False
    return True
