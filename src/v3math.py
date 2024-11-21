# Here are some basic operations for vectors

def v3dpp(x1, x2):
    # distance between a point and a point
    out = 0.0
    for ix in range(3):
        out += (x1[ix]-x2[ix])*(x1[ix]-x2[ix])
    return out**0.5


def v3matchpp(x1, x2, tol=1e-6):
    # if two points are the same
    for ix in range(3):
        if abs(x1[ix]-x2[ix]) > tol:
            return False
    return True


def v3pvm(x, M):
    # product of a vector and a matrix
    # like apc = apd * lc
    xout = [0.0, 0.0, 0.0]
    for ix1 in range(3):
        for ix2 in range(3):
            xout[ix2] += x[ix1]*M[ix1][ix2]
    return xout
