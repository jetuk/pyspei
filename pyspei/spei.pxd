

cdef extern from "spei.c":
    #void spei(float dataSeries[], int n, int seasons, float speiSeries[])
    void thornthwaite(float tempSeries[], int n, float lat, float etpSeries[])
    void upward(float series[], int n)
    void pwm(float series[], int n, float beta[], float A, float B, int isBeta)
    void logLogisticFit(float beta[], float logLogisticParams[])
    float logLogisticCDF(float value, float params[])
    float standardGaussianInvCDF(float prob)