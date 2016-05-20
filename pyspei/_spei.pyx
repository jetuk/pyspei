
cimport spei as c_spei
import numpy as np
cimport numpy as np

cpdef spei(float[:] balanceSeries, int months, int seasons, float[:] speiSeries):
    cdef int n = balanceSeries.shape[0]
    cdef int i, j

	# Compute the cumulative series
    cdef int acumN = n-months+1
    cdef float[:] acumSeries = np.zeros(acumN, dtype=np.float32)
    for i in range(months-1, n):
        for j in range(months):
            acumSeries[i-months+1] += balanceSeries[i-j]
    # Calculate SPEI
    c_spei.spei(&acumSeries[0], acumN, seasons, &speiSeries[0])

cpdef thornthwaite(float[:] tempSeries, float lat, float[:] etpSeries):

    cdef int n = tempSeries.shape[0]
    c_spei.thornthwaite(&tempSeries[0], n, lat, &etpSeries[0])

