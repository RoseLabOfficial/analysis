import numpy as np

LOG2PI = np.log(2.0 * np.pi)


def _as_pvec(x, P, name):
    """Ensure x is shape (P,) or scalar -> broadcast to (P,)."""
    x = np.asarray(x, dtype=float)
    if x.ndim == 0:
        return np.full(P, float(x))
    if x.shape != (P,):
        raise ValueError(f"{name} must be scalar or shape (P,), got {x.shape}")
    return x


def posterior_C_batch_grid(
    a, b, C_grid,
    sigma_A, sigma_B,
    mu_C, tau_C,
    mu_B, tau_B,
):
    """
    Vectorized posterior p(C | a,b) for P independent problems.

    Inputs
    ------
    a, b : arrays shape (P, N)
        Observations (can contain np.nan for missing pairs).
    C_grid : array shape (M,)
        Common grid of C values to evaluate.
    sigma_A, sigma_B : scalar or (P,)
    mu_C, tau_C : scalar or (P,)
    mu_B, tau_B : scalar or (P,)

    Returns
    -------
    pdf : array shape (P, M)
        Normalized posterior density over C_grid for each problem.
    logpost : array shape (P, M)
        Unnormalized log posterior (up to additive constant) for each problem.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape or a.ndim != 2:
        raise ValueError("a and b must be 2D arrays with the same shape (P, N).")
    P, N = a.shape

    C = np.asarray(C_grid, dtype=float)
    if C.ndim != 1:
        raise ValueError("C_grid must be 1D.")
    M = C.size

    # Broadcast parameters to (P,)
    sigma_A = _as_pvec(sigma_A, P, "sigma_A")
    sigma_B = _as_pvec(sigma_B, P, "sigma_B")
    mu_C    = _as_pvec(mu_C,    P, "mu_C")
    tau_C   = _as_pvec(tau_C,   P, "tau_C")
    mu_B    = _as_pvec(mu_B,    P, "mu_B")
    tau_B   = _as_pvec(tau_B,   P, "tau_B")

    # Mask: only use rows where BOTH a and b are present
    mask = np.isfinite(a) & np.isfinite(b)
    n = mask.sum(axis=1).astype(float)  # shape (P,)

    # Handle pathological case: some problems have zero valid observations
    # (Posterior will just equal the prior on that row.)
    # We'll still compute safely by setting stats to 0 when n=0.
    a0 = np.where(mask, a, 0.0)
    b0 = np.where(mask, b, 0.0)

    Sa  = a0.sum(axis=1)           # sum a
    Sb  = b0.sum(axis=1)           # sum b
    Sa2 = (a0 * a0).sum(axis=1)    # sum a^2
    Sb2 = (b0 * b0).sum(axis=1)    # sum b^2
    Sab = (a0 * b0).sum(axis=1)    # sum a*b

    # Shapes for broadcasting
    # params: (P,1), C: (1,M)
    C_row = C[None, :]                 # (1,M)
    muB   = mu_B[:, None]              # (P,1)
    varB  = (tau_B ** 2)[:, None]      # (P,1)
    sigA2 = (sigma_A ** 2)[:, None]    # (P,1)
    sigB2 = (sigma_B ** 2)[:, None]    # (P,1)
    ncol  = n[:, None]                 # (P,1)

    # Covariance components for [a,b] | C with B integrated out
    cov11 = (C_row**2) * varB + sigA2          # (P,M)
    cov12 = C_row * varB                       # (P,M)
    cov22 = varB + sigB2                       # (P,1) broadcastable
    cov22 = np.broadcast_to(cov22, (P, M))     # (P,M)

    det = cov11 * cov22 - cov12 * cov12        # (P,M)
    # Numerical safety: det should be >0; clamp tiny negatives from roundoff
    det = np.maximum(det, 1e-300)
    logdet = np.log(det)

    # Sufficient-statistics forms for residual sums
    # ra_i = a_i - C*muB, rb_i = b_i - muB
    # Sum ra^2 depends on C; Sum rb^2 constant; Sum ra*rb depends on C.
    Sa  = Sa[:, None]
    Sb  = Sb[:, None]
    Sa2 = Sa2[:, None]
    Sb2 = Sb2[:, None]
    Sab = Sab[:, None]

    Srb2 = Sb2 - 2.0 * muB * Sb + ncol * (muB ** 2)  # (P,1)

    Sra2 = Sa2 - 2.0 * C_row * muB * Sa + ncol * (C_row ** 2) * (muB ** 2)  # (P,M)

    Srarb = (
        Sab
        - muB * Sa
        - C_row * muB * Sb
        + ncol * C_row * (muB ** 2)
    )  # (P,M)

    # Quadratic form summed over i: (cov22*Sum ra^2 - 2cov12*Sum ra rb + cov11*Sum rb^2)/det
    quad_sum = (cov22 * Sra2 - 2.0 * cov12 * Srarb + cov11 * Srb2) / det  # (P,M)

    # Log-likelihood: sum_i log N_2( [a_i,b_i] | mean(C), Sigma(C) )
    # = -0.5 * [ n*(2log2pi + logdet) + quad_sum ]
    ll = -0.5 * (ncol * (2.0 * LOG2PI + logdet) + quad_sum)  # (P,M)

    # Prior on C: N(mu_C, tau_C^2)
    muC  = mu_C[:, None]
    varC = (tau_C ** 2)[:, None]
    lpC = -0.5 * (np.log(2.0 * np.pi * varC) + ((C_row - muC) ** 2) / varC)  # (P,M)

    logpost = lpC + ll

    # If n=0 for a row, ll should be 0 (no data); current formula yields ll=0 because n=0
    # BUT quad_sum might have 0/ det artifacts because stats are 0; it's still 0.
    # To be explicit, we can force ll=0 when n=0:
    zero_data = (n == 0.0)
    if np.any(zero_data):
        logpost[zero_data, :] = lpC[zero_data, :]

    # Normalize each row to a density over C_grid via trapezoid integral
    maxlp = np.max(logpost, axis=1)
    w = np.exp(logpost - maxlp[:, np.newaxis])
    Z = np.trapezoid(w, C, axis=1)
    pdf = w / Z[:, np.newaxis]

    return pdf, logpost


def summarize_posterior_batch(C_grid, pdf, cred_mass=0.95):
    """
    Vectorized posterior summaries for pdf shape (P,M) over C_grid shape (M,).
    Returns dict of arrays shape (P,).
    """
    C = np.asarray(C_grid, dtype=float)
    pdf = np.asarray(pdf, dtype=float)
    if pdf.ndim != 2 or C.ndim != 1 or pdf.shape[1] != C.size:
        raise ValueError("pdf must be (P,M) and C_grid must be (M,)")

    # Mean/var via integration
    mean = np.trapezoid(pdf * C[None, :], C, axis=1)
    second = np.trapezoid(pdf * (C[None, :] ** 2), C, axis=1)
    var = np.maximum(0.0, second - mean**2)
    std = np.sqrt(var)

    # CDF via cumulative trapezoids
    dx = np.diff(C)
    cdf = np.zeros_like(pdf)
    cdf[:, 1:] = np.cumsum(0.5 * (pdf[:, 1:] + pdf[:, :-1]) * dx[None, :], axis=1)

    def q(p):
        # row-wise interpolation
        return np.array([np.interp(p, cdf_i, C) for cdf_i in cdf])

    alpha = (1.0 - cred_mass) / 2.0
    lo = q(alpha)
    med = q(0.5)
    hi = q(1.0 - alpha)

    map_idx = np.argmax(pdf, axis=1)
    mapC = C[map_idx]

    return {
        "mean": mean,
        "std": std,
        "median": med,
        "map": mapC,
        f"ci_{int(cred_mass*100)}_lo": lo,
        f"ci_{int(cred_mass*100)}_hi": hi,
    }


if __name__ == "__main__":
    # Example: P=3 problems, N=5 obs each
    a = np.array([
        [10.2,  9.7, 10.5,  np.nan, 10.1],
        [ 3.1,  2.9,  3.3,   3.0,    3.2],
        [20.0, 19.5, 20.4,  20.2,   19.8],
    ])
    b = np.array([
        [2.1, 1.9, 2.2, np.nan, 2.0],
        [0.7, 0.8, 0.75, 0.72, 0.78],
        [4.0, 3.9, 4.1, 4.05, 3.95],
    ])

    C_grid = np.linspace(-5, 15, 4001)

    pdf, logpost = posterior_C_batch_grid(
        a, b, C_grid,
        sigma_A=0.5, sigma_B=0.2,
        mu_C=5.0, tau_C=2.0,
        mu_B=np.array([2.0, 0.75, 4.0]),  # per-problem B prior mean
        tau_B=1.0
    )

    summ = summarize_posterior_batch(C_grid, pdf, cred_mass=0.95)
    for k, v in summ.items():
        print(k, v)
