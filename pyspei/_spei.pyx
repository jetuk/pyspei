
cimport spei as c_spei
import numpy as np
cimport numpy as np

cpdef spei(float[:] balanceSeries, int months, int seasons, float[:] speiSeries):
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
    cdef int n = balanceSeries.shape[0]
    cdef int i, j

	# Compute the cumulative series
    cdef int acumN = n-months+1
    cdef float[:] acumSeries = np.zeros(acumN, dtype=np.float32)
    for i in range(months-1, n):
        for j in range(months):
            acumSeries[i-months+1] += balanceSeries[i-j]
    # Calculate SPEI
    #c_spei.spei(&acumSeries[0], acumN, seasons, &speiSeries[0])
    _spei(acumSeries, acumN, seasons, speiSeries)


cpdef thornthwaite(float[:] tempSeries, float lat, float[:] etpSeries):

    cdef int n = tempSeries.shape[0]
    c_spei.thornthwaite(&tempSeries[0], n, lat, &etpSeries[0])

cdef _spei(float[:] dataSeries, int n, int seasons, float[:] speiSeries):
    """
    Python port of C spei() function.

    Calculates the Standardized Precipitation-Evapotransporation Index
     from a series of climatic balance (precipitation minus etp). The
     SPEI is the standardized value of the climatic balance (P-ETP),
     computed following a Log Logistic probability distribution.

    """
    cdef int i, j, k
    cdef int nSeason = n // seasons + 1
    cdef float[:] seasonSeries = np.empty(nSeason, dtype=np.float32)
    cdef float[:] beta = np.empty(3, dtype=np.float32)
    # NB the 1 based indexing.
    cdef float[:, :] logLogisticParams = np.empty((seasons+1, 3), dtype=np.float32)

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

        # Calculate the standardized values
        for i in range(j-1, n, seasons):
            speiSeries[i] = c_spei.logLogisticCDF(dataSeries[i], &logLogisticParams[j, 0])
            speiSeries[i] = -c_spei.standardGaussianInvCDF(speiSeries[i])
