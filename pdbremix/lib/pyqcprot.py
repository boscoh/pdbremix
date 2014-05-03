import math


def make_correlation_matrix(coords1, coords2):
    N = len(coords1)
    assert N == len(coords2)

    G1 = 0.0
    G2 = 0.0
    
    A = [0.0 for i in range(9)]
    
    for i in xrange(N):
        x1 = coords1[i][0]
        y1 = coords1[i][1]
        z1 = coords1[i][2]

        G1 += (x1*x1 + y1*y1 + z1*z1)

        x2 = coords2[i][0]
        y2 = coords2[i][1]
        z2 = coords2[i][2]

        G2 += (x2*x2 + y2*y2 + z2*z2)

        A[0] +=  (x1 * x2)
        A[1] +=  (x1 * y2)
        A[2] +=  (x1 * z2)

        A[3] +=  (y1 * x2)
        A[4] +=  (y1 * y2)
        A[5] +=  (y1 * z2)

        A[6] +=  (z1 * x2)
        A[7] +=  (z1 * y2)
        A[8] +=  (z1 * z2)
        
    E0 = (G1 + G2) * 0.5

    return E0, A


def calc_rms_rot(coords1, coords2):
    E0, A = make_correlation_matrix(coords1, coords2)
    N = len(coords1)

    oldg = 0.0
    evecprec = 1e-6
    evalprec = 1e-14
    
    Sxx = A[0] 
    Sxy = A[1]
    Sxz = A[2]
    Syx = A[3]
    Syy = A[4]
    Syz = A[5]
    Szx = A[6]
    Szy = A[7]
    Szz = A[8]

    Sxx2 = Sxx * Sxx
    Syy2 = Syy * Syy
    Szz2 = Szz * Szz

    Sxy2 = Sxy * Sxy
    Syz2 = Syz * Syz
    Sxz2 = Sxz * Sxz

    Syx2 = Syx * Syx
    Szy2 = Szy * Szy
    Szx2 = Szx * Szx

    SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz)
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2

    C = [0.0 for i in range(3)]

    C[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2)
    C[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz)

    SxzpSzx = Sxz + Szx
    SyzpSzy = Syz + Szy
    SxypSyx = Sxy + Syx
    SyzmSzy = Syz - Szy
    SxzmSzx = Sxz - Szx
    SxymSyx = Sxy - Syx
    SxxpSyy = Sxx + Syy
    SxxmSyy = Sxx - Syy
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2

    C[0] = (Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
         + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
         + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
         + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
         + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
         + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz)))

    mxEigenV = E0
    n_iter = 50
    for i in range(n_iter):
        oldg = mxEigenV
        x2 = mxEigenV*mxEigenV
        b = (x2 + C[2])*mxEigenV
        a = b + C[1]
        delta = ((a*mxEigenV + C[0])/(2.0*x2*mxEigenV + b + a))
        mxEigenV -= delta
        if (abs(mxEigenV - oldg) < abs((evalprec)*mxEigenV)):
            break
    else:
        raise Exception("More iterations needed to find eigenvalue")

    val = 2.0 * (E0 - mxEigenV)/float(N)
    if abs(val) < evecprec:
        rms = 0.0
    else:
        rms = math.sqrt(val)

    rot = [0.0 for i in range(9)]

    a11 = SxxpSyy + Szz-mxEigenV
    a12 = SyzmSzy
    a13 = - SxzmSzx
    a14 = SxymSyx

    a21 = SyzmSzy
    a22 = SxxmSyy - Szz-mxEigenV
    a23 = SxypSyx
    a24= SxzpSzx

    a31 = a13
    a32 = a23
    a33 = Syy-Sxx-Szz - mxEigenV
    a34 = SyzpSzy

    a41 = a14
    a42 = a24
    a43 = a34
    a44 = Szz - SxxpSyy - mxEigenV

    a3344_4334 = a33 * a44 - a43 * a34
    a3244_4234 = a32 * a44-a42*a34
    a3243_4233 = a32 * a43 - a42 * a33
    a3143_4133 = a31 * a43-a41*a33
    a3144_4134 = a31 * a44 - a41 * a34
    a3142_4132 = a31 * a42-a41*a32

    q1 =  a22*a3344_4334-a23*a3244_4234+a24*a3243_4233
    q2 = -a21*a3344_4334+a23*a3144_4134-a24*a3143_4133
    q3 =  a21*a3244_4234-a22*a3144_4134+a24*a3142_4132
    q4 = -a21*a3243_4233+a22*a3143_4133-a23*a3142_4132

    qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4

    # The following code tries to calculate another column in the adjoint
    # matrix when the norm of the current column is too small. Usually
    # this commented block will never be activated.  To be absolutely safe
    # this should be uncommented, but it is most likely unnecessary.

    if (qsqr < evecprec):
        q1 =  a12*a3344_4334 - a13*a3244_4234 + a14*a3243_4233
        q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133
        q3 =  a11*a3244_4234 - a12*a3144_4134 + a14*a3142_4132
        q4 = -a11*a3243_4233 + a12*a3143_4133 - a13*a3142_4132
        qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4

        if (qsqr < evecprec):
            a1324_1423 = a13 * a24 - a14 * a23 
            a1224_1422 = a12 * a24 - a14 * a22
            a1223_1322 = a12 * a23 - a13 * a22
            a1124_1421 = a11 * a24 - a14 * a21
            a1123_1321 = a11 * a23 - a13 * a21
            a1122_1221 = a11 * a22 - a12 * a21

            q1 =  a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322
            q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321
            q3 =  a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221
            q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221
            qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4

            if (qsqr < evecprec):
                q1 =  a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322
                q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321
                q3 =  a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221
                q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221
                qsqr = q1*q1 + q2 *q2 + q3*q3 + q4*q4

                if (qsqr < evecprec):
                    # if qsqr is still too small, return the identity matrix. #
                    rot[0] = rot[4] = rot[8] = 1.0
                    rot[1] = rot[2] = rot[3] = rot[5] = rot[6] = rot[7] = 0.0

                    return rms, rot


    normq = math.sqrt(qsqr)
    q1 /= normq
    q2 /= normq
    q3 /= normq
    q4 /= normq

    a2 = q1 * q1
    x2 = q2 * q2
    y2 = q3 * q3
    z2 = q4 * q4

    xy = q2 * q3
    az = q1 * q4
    zx = q4 * q2
    ay = q1 * q3
    yz = q3 * q4
    ax = q1 * q2

    rot[0] = a2 + x2 - y2 - z2
    rot[1] = 2 * (xy + az)
    rot[2] = 2 * (zx - ay)
    rot[3] = 2 * (xy - az)
    rot[4] = a2 - x2 + y2 - z2
    rot[5] = 2 * (yz + ax)
    rot[6] = 2 * (zx + ay)
    rot[7] = 2 * (yz - ax)
    rot[8] = a2 - x2 - y2 + z2

    return rms, rot


      
