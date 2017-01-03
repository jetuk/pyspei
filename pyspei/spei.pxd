

cdef extern from "spei.c":
    #void spei(double dataSeries[], int n, int seasons, double speiSeries[])
    void thornthwaite(double tempSeries[], int n, double lat, double etpSeries[])
    void upward(double series[], int n)
    void pwm(double series[], int n, double beta[], double A, double B, int isBeta)
    void logLogisticFit(double beta[], double logLogisticParams[])
    double logLogisticCDF(double value, double params[])
    double standardGaussianInvCDF(double prob)