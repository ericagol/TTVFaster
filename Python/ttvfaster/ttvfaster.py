# This routine computes the first-order
# transit timing variations in Agol & Deck (2015).  Please
# cite the paper if you make use of this code in your research.

import numpy as np
import matplotlib.pyplot as plt


class Planet(object):
    def __init__(self, mass_ratio=None, trans0=None, period=None, ecos=None,
                 esin=None):
        self.mass_ratio = mass_ratio
        self.trans0 = trans0
        self.period = period
        self.ecos = ecos
        self.esin = esin


def call_ttv(jmax):
    """
    This routine gives an example of a call of compute_ttv.py
    which computes the first-order eccentricity TTVs, from
    Agol & Deck (2015).
    It uses parameters appropriate for the outer two planets
    of Kepler-62e/f.

    Parameters
    ----------
    jmax:  maximum j to evaluate
    """
    data = np.loadtxt('kepler62ef_planets.txt', delimiter=',')
    # Compute 40 transits of inner and 20 of outer planet
    nt1 = 40
    nt2 = 20
    # Set up uniform ephemeris to compute TTVs
    time1 = data[2] + np.arange(nt1) * data[1]
    time2 = data[7] + np.arange(nt2) * data[6]
    # Rearrange planet parameters for ingest to compute_ttv.py
    param = np.array(data)
    param[1] = data[2]
    param[2] = data[1]
    param[6] = data[7]
    param[7] = data[6]

    # Model contains the transit times
    model1, model2 = compute_ttv(time1, time2, param, jmax)

    # Subtract the transit times to only plot the TTVs
    plt.plot(time1, model1 - time1)
    plt.plot(time2, model2 - time2)
    plt.show()
    np.savetxt('inner_ttv.txt', model1)
    np.savetxt('outer_ttv.txt', model2)


def compute_ttv(time1, time2, param, jmax):
    """
    Computes transit-timing variations to linear order in
    eccentricity for non-resonant, plane-parallel planets
    to first order in mass ratio and eccentricity.

    Parameters
    ----------
    time1:  contains the times of transit of inner planet
    time2:  contains the times of transit of outer planet
    jmax:  maximum j to evaluate
    param:  contains the (mass ratio,t0,period,e*cos(omega),
    e*sin(omega)) for each planet

    Returns
    -------
    model1:  Contains the transit timing model with the first
    set for the inner planet.
    model2:  Contains the transit timing model with the first
    set for the outer planet.
    """
    # Insert the parameters of the model into a structure for each planet
    p1 = Planet(mass_ratio=param[0], trans0=param[1], period=param[2],
                ecos=param[3], esin=param[4])
    p2 = Planet(mass_ratio=param[5], trans0=param[6], period=param[7],
                ecos=param[8], esin=param[9])

    twopi = 2.0 * np.pi
    # Compute the semi-major axis ratio of the planets
    alpha = abs(p1.period / p2.period) ** (2. / 3.)
    # Compute the longitudes of the planets at times of transit of
    # planet 1 (equation 32)
    lam11 = twopi * (time1 - p1.trans0) / p1.period + 2 * p1.esin
    lam21 = twopi * (time1 - p2.trans0) / p2.period + 2 * p2.esin
    # Compute the longitudes of the planets at times of transit of
    # planet 2 (equation 32)
    lam12 = twopi * (time2 - p1.trans0) / p1.period + 2 * p1.esin
    lam22 = twopi * (time2 - p2.trans0) / p2.period + 2 * p2.esin
    # Compute difference in longitudes at times of transit of planets 1 and 2
    psi1 = lam11 - lam21
    psi2 = lam12 - lam22
    # Compute the coefficients (need one higher than jmax)
    f1, f2 = ttv_succinct(jmax + 1, alpha)

    ttv1 = np.zeros(time1.size)
    ttv2 = np.zeros(time2.size)
    # Compute TTVs for inner planet (equation 33)
    for j in range(1, jmax + 1):
        ttv1 += f1[j, 0] * np.sin(j * psi1)
        ttv1 += f1[j, 1] * (p1.ecos * np.sin((j - 1) * lam11 - j * lam21) +
                            p1.esin * np.cos((j - 1) * lam11 - j * lam21))
        ttv1 += f1[j, 2] * (p1.ecos * np.sin((j + 1) * lam11 - j * lam21) -
                            p1.esin * np.cos((j + 1) * lam11 - j * lam21))
        ttv1 += f1[j - 1, 3] * (p2.ecos * np.sin((j - 1) * lam11 - j * lam21) +
                                p2.esin * np.cos((j - 1) * lam11 - j * lam21))
        ttv1 += f1[j + 1, 4] * (p2.ecos * np.sin((j + 1) * lam11 - j * lam21) -
                                p2.esin * np.cos((j + 1) * lam11 - j * lam21))
    # Now multiply by the mass ratio and divide by mean motion
    ttv1 = ttv1 * p1.period * p2.mass_ratio / twopi

    # Compute TTVs for outer planet (equation 33)
    for j in range(1, jmax + 1):
        ttv2 += f2[j, 0] * np.sin(j * psi2)
        ttv2 += f2[j, 1] * (p2.ecos * np.sin(j * lam12 - (j + 1) * lam22) +
                            p2.esin * np.cos(j * lam12 - (j + 1) * lam22))
        ttv2 += f2[j, 2] * (p2.ecos * np.sin(j * lam12 - (j - 1) * lam22) -
                            p2.esin * np.cos(j * lam12 - (j - 1) * lam22))
        ttv2 += f2[j + 1, 3] * (p1.ecos * np.sin(j * lam12 - (j + 1) * lam22) +
                                p1.esin * np.cos(j * lam12 - (j + 1) * lam22))
        ttv2 += f2[j - 1, 4] * (p1.ecos * np.sin(j * lam12 - (j - 1) * lam22) -
                                p1.esin * np.cos(j * lam12 - (j - 1) * lam22))
    # Now multiply by the mass ratio and divide by mean motion
    ttv2 = ttv2 * p2.period * p1.mass_ratio / twopi

    # Add the TTVs to the ephemeris and return to user
    model1 = p1.trans0 + np.arange(time1.size) * p1.period + ttv1
    model2 = p2.trans0 + np.arange(time2.size) * p2.period + ttv2
    # The following lines can be used if only the TTVs are desired
    #    model1= ttv1
    #    model2= ttv2
    return model1, model2


# Define u & v functions (equation 34)
def u(gamma, c1, c2):
    return (((3.0 + gamma ** 2) * c1 + 2.0 * gamma * c2) /
            gamma ** 2 / (1.0 - gamma ** 2))


def v(zeta, d1, d2, m):
    # m = +/-1
    return ((m * (1.0 - zeta ** 2) + 6.0 * zeta) * d1 +
            (2.0 + zeta ** 2) * d2) / (zeta * (1.0 - zeta ** 2) *
                                       (zeta + m) * (zeta + 2.0 * m))


# The following routine computes the coefficients for the first-order
# transit timing variations in Agol & Deck (2015).  Please
# cite the paper if you make use of this code in your research.

def ttv_succinct(jmax, alpha):
    """
    Succinct form for coefficients of first-order TTV formula

    Parameters
    ----------
    jmax:  maximum value of j over which to compute the coefficients
    alpha:  a_1/a_2 of the two planets

    Returns
    -------
    b:  Matrix of Laplace coefficients.
    f1:  Coefficients for the inner planet.  For each value of
        j=0 to jmax, there are 5 coefficients:
        the f_{1,j}^(0), f_{1,j}^(+-1), and
        f_{1,j}^(+-2) coefficients.

        The +-1 coefficients correspond to
        j*(lam_1-lam_2)+-(lam_1-omega_1) arguments.

        The +-2 coefficients correspond to
        j*(lam_1-lam_2)+-(lam_1-omega_2) arguments.
    f2:  Coefficients for the outer planet.  For each value of
        j=0 to jmax, there are 5 coefficients:
        the f_{1,j}^(0), f_{1,j}^(+-2), and
        f_{1,j}^(+-1) coefficients.

        The +-1 coefficients correspond to
        j*(lam_1-lam_2)+-(lam_2-omega_1) arguments.

        The +-2 coefficients correspond to
        j*(lam_1-lam_2)+-(lam_2-omega_2) arguments.
    """
    f1 = np.zeros((jmax + 1, 5))
    f2 = np.zeros((jmax + 1, 5))

    # Compute the Laplace coefficients
    b = laplace_coefficients3(jmax, alpha)
    # Now loop over j values
    for j in range(jmax + 1):
        if j == 1:
            dj1 = 1.0
        else:
            dj1 = 0.0
        beta = j * (1 - alpha ** 1.5)
        kappa = beta / alpha ** 1.5

        # Compute the disturbing function coefficients A_jmn (equation 31)
        A_j00 = b[j, 0]
        A_j10 = alpha * b[j, 1]
        A_j01 = -(A_j10 + A_j00)
        A_j20 = alpha ** 2 * b[j, 2]
        A_j11 = -(2 * A_j10 + A_j20)
        A_j02 = 2 * A_j00 + 4 * A_j10 + A_j20

        # Inner planet coefficients, in order k=0,-1,1,-2,2 (see Table 1)
        gamma = beta + np.array([0, -1, 1, -alpha ** 1.5, alpha ** 1.5])

        c1 = alpha * j * (A_j00 * np.array([1, -j, j, j, -j]) -
                          .5 * A_j01 * np.array([0, 0, 0, 1, 1]) -
                          .5 * A_j10 * np.array([0, 1, 1, 0, 0]) -
                          alpha * dj1 * np.array([1, -1.5, .5, 2, 0]))
        c2 = alpha * (A_j10 * np.array([1, -j, j, j, -j]) -
                      .5 * A_j11 * np.array([0, 0, 0, 1, 1]) -
                      .5 * A_j20 * np.array([0, 1, 1, 0, 0]) -
                      alpha * dj1 * np.array([1, -1, 1, 2, 0]))
        if j >= 2:
            for k in range(5):
                f1[j, k] = u(gamma[k], c1[k], c2[k])
        else:
            if j == 0:
                f1[j, 3] = u(gamma[3], c1[3], c2[3])
            else:
                for k in range(4):
                    f1[j, k] = u(gamma[k], c1[k], c2[k])

        # Add in the k=\pm 1 coefficients (note that d1 & d2 are the same as
        # c1 & c2 for k=0)
        if j >= 1:
            f1[j, 1] += v(beta, c1[0], c2[0], -1)
            f1[j, 2] += v(beta, c1[0], c2[0], 1)

        # Now for the outer planet
        # Outer planet coefficients, in order k=0,-2,2,-1,1 (see Table 1)
        gamma = kappa + np.array([0, -1, 1, -1.0 / alpha ** 1.5,
                                  1.0 / alpha ** 1.5])
        c1 = -j * (A_j00 * np.array([1, j, -j, -j, j]) -
                   .5 * A_j01 * np.array([0, 1, 1, 0, 0]) -
                   .5 * A_j10 * np.array([0, 0, 0, 1, 1]) -
                   dj1 / alpha ** 2 * np.array([1, .5, -1.5, 0, 2]))
        c2 = (A_j01 * np.array([1, j, -j, -j, j]) -
              .5 * A_j11 * np.array([0, 0, 0, 1, 1]) -
              .5 * A_j02 * np.array([0, 1, 1, 0, 0]) -
              dj1 / alpha ** 2 * np.array([1, 1, -1, 0, 2]))
        if j >= 2:
            for k in range(4):
                f2[j, k] = u(gamma[k], c1[k], c2[k])
        else:
            if j == 1:
                for k in range(3):
                    f2[j, k] = u(gamma[k], c1[k], c2[k])
        f2[j, 4] = u(gamma[4], c1[4], c2[4])
        # Add in the k=\pm 2 coefficients (note that d1 & d2 are the same as
        # c1 & c2 for k=0)
        if j >= 1:
            f2[j, 1] += v(kappa, c1[0], c2[0], -1)
            f2[j, 2] += v(kappa, c1[0], c2[0], 1)

    return f1, f2


def laplace_coefficients3(j, alpha):
    # This computes the Laplace coefficients via recursion.
    # Compute the highest two Laplace coefficients using
    # Wisdom's series approach
    b = np.zeros((j + 1, 3))
    b[j, 0] = laplace_wisdom(.5, j, 0, alpha)
    b[j - 1, 0] = laplace_wisdom(.5, j - 1, 0, alpha)
    b[j, 1] = laplace_wisdom(.5, j, 1, alpha) / alpha
    b[j - 1, 1] = laplace_wisdom(.5, j - 1, 1, alpha) / alpha
    b[j, 2] = laplace_wisdom(.5, j, 2, alpha) / alpha / alpha
    b[j - 1, 2] = laplace_wisdom(.5, j - 1, 2, alpha) / alpha / alpha

    # The rest can be found with the recursion formulae
    j0 = j - 2
    while j0 >= 0:
        # Recurrence relations (derived on 9/9/14; also see Brouwer & Clemence)
        jd = float(j0)
        b[j0, 0] = (-(2.0 * jd + 3.0) * b[j0 + 2, 0] +
                    (1.0 + alpha * alpha) / alpha * 2.0 * (jd + 1.0) *
                    b[j0 + 1, 0]) / (2.0 * jd + 1.0)
        # Recurrence for first derivative
        b[j0, 1] = b[j0 + 2, 1] - (2.0 * alpha * (jd + 1.0) * b[j0 + 1, 0] -
                                   jd * b[j0, 0] -
                                   (jd + 2.0) * b[j0 + 2, 0]) / alpha
        # Recurrence for second derivative
        b[j0, 2] = b[j0 + 2, 2] - 2.0 * (jd + 1.0) * b[j0 + 1, 1] - \
                   (jd * b[j0, 0] +
                    (jd + 2.0) * b[j0 + 2, 0]) / alpha / alpha + \
                   (jd * b[j0, 1] + (jd + 2.0) * b[j0 + 2, 1]) / alpha
        j0 -= 1
    return b


def laplace_wisdom(s, i, j, a):
    """
    Code due to Jack Wisdom.
    Compute Laplace coefficients and Leverrier derivatives
          j
     j   d     i
    a   ---   b (a)
          j    s
        da

    by series summation.
    """
    LAPLACE_EPS = 1.0e-12

    asquare = a * a

    i = abs(i)

    # compute first term in sum
    if j <= i:
        factor4 = 1.0
        for k in range(j):
            factor4 *= float(i - k)
        isum = factor4
        q0 = 0
    else:
        q0 = (j + 1 - i) / 2
        isum = 0.
        factor4 = 1.

    # compute factors for terms in sum
    factor1 = s
    factor2 = s + float(i)
    factor3 = float(i) + 1.0
    q = 1
    # no contribution for q = 0
    for q in range(1, q0):
        factor1 *= s + float(q)
        factor2 *= s + float(i + q)
        factor3 *= float(i + q + 1)

    term = asquare * factor1 * factor2 / (factor3 * float(q))

    #  sum series
    while (term * factor4) > LAPLACE_EPS:
        factor4 = 1.0
        for k in range(j):
            factor4 *= 2.0 * float(q) + float(i - k)
        isum += term * factor4
        factor1 += 1.0
        factor2 += 1.0
        factor3 += 1.0
        q += 1
        term = term * asquare * factor1 * factor2 / (factor3 * float(q))

    # fix coefficient
    for k in range(i):
        isum = isum * (s + float(k)) / (float(k) + 1.0)

    if q0 <= 0:
        isum = isum * 2.0 * a ** i
    else:
        isum = isum * 2.0 * a ** (2 * q0 + i - 2)

    return isum
