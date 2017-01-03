
cimport spei as c_spei
import numpy as np
cimport numpy as np


cdef double[:] acummulate_series(double[:] balanceSeries, int months):
    cdef int n = balanceSeries.shape[0]
    cdef int i, j
	# Compute the cumulative series
    cdef int acumN = n-months+1
    cdef double[:] acumSeries = np.zeros(acumN, dtype=np.float64)
    for i in range(months-1, n):
        for j in range(months):
            acumSeries[i-months+1] += balanceSeries[i-j]

    return acumSeries

cpdef spei(double[:] balanceSeries, int months, int seasons, double[:] speiSeries,
           double[:] beta=None, double[:, :] logLogisticParams=None):
    """
    Calculate SPEI for a monthly series of effective rainfall

    Parameters
    ----------
    balanceSeries
        Monthly effective raninfall array
    months
        Size of the window in months
    seasons
        Should be 12 according to SPEI documentation.
    speiSeries
        Preallocated output array which is populated with SPEI values

    Returns
    -------
    None

    """
    cdef double[:] acumSeries = acummulate_series(balanceSeries, months)
    cdef int acumN = acumSeries.shape[0]

    if beta is None and logLogisticParams is None:
        beta = np.empty(3, dtype=np.float64)
        # NB the 1 based indexing.
        logLogisticParams = np.empty((seasons+1, 3), dtype=np.float64)
        # Calculate log-logistic distribution parameters
        _fit_loglogistic(acumSeries, seasons, beta, logLogisticParams)

    # Calculate SPEI
    _spei(acumSeries,  seasons,  beta, logLogisticParams, speiSeries)
    return beta, logLogisticParams


def fit_loglogistic(double[:] balanceSeries, int months, int seasons):
    cdef double[:] acumSeries = acummulate_series(balanceSeries, months)
    cdef int acumN = acumSeries.shape[0]

    cdef double[:] beta = np.empty(3, dtype=np.float64)
    # NB the 1 based indexing.
    cdef double[:, :] logLogisticParams = np.zeros((seasons+1, 3), dtype=np.float64)

    # Calculate log-logistic distribution parameters
    _fit_loglogistic(acumSeries, seasons, beta, logLogisticParams)
    return beta, logLogisticParams

cpdef thornthwaite(double[:] tempSeries, double lat, double[:] etpSeries):

    cdef int n = tempSeries.shape[0]
    c_spei.thornthwaite(&tempSeries[0], n, lat, &etpSeries[0])


cdef _fit_loglogistic(double[:] dataSeries, int seasons, double[:] beta, double[:, :] logLogisticParams):
    """
    Python port of first half of C spei() function.

    Calculates the Standardized Precipitation-Evapotransporation Index
     from a series of climatic balance (precipitation minus etp). The
     SPEI is the standardized value of the climatic balance (P-ETP),
     computed following a Log Logistic probability distribution.

    """
    cdef int i, j, k
    cdef int n = dataSeries.shape[0]
    cdef int nSeason = n // seasons + 1
    cdef double[:] seasonSeries = np.empty(nSeason, dtype=np.float64)

	# Loop through all seasons defined by seasons
    for j in range(1, seasons+1):
		# Extract and sort the seasonal series
        k = 0
        for i in range(j-1, n, seasons):
            seasonSeries[k] = dataSeries[i]
            k+=1

        nSeason = k
        c_spei.upward(&seasonSeries[0], nSeason)
        # Compute probability weighted moments
        c_spei.pwm(&seasonSeries[0], nSeason, &beta[0], 0, 0, 0)
        # Fit a Log Logistic probability function
        c_spei.logLogisticFit(&beta[0], &logLogisticParams[j, 0])


cdef _spei(double[:] dataSeries, int seasons, double[:] beta, double[:, :] logLogisticParams, double[:] speiSeries):
    """
    Python port of second half of C spei() function.

    Calculates the Standardized Precipitation-Evapotransporation Index
     from a series of climatic balance (precipitation minus etp). The
     SPEI is the standardized value of the climatic balance (P-ETP),
     computed following a Log Logistic probability distribution.

    """
    cdef int n = dataSeries.shape[0]
    cdef int i, j
	# Loop through all seasons defined by seasons
    for j in range(1, seasons+1):
        # Calculate the standardized values
        for i in range(j-1, n, seasons):
            speiSeries[i] = c_spei.logLogisticCDF(dataSeries[i], &logLogisticParams[j, 0])
            speiSeries[i] = -c_spei.standardGaussianInvCDF(speiSeries[i])
