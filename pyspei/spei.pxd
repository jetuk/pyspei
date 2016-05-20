

cdef extern from "spei.c":
    void spei(float dataSeries[], int n, int seasons, float speiSeries[])
    void thornthwaite(float tempSeries[], int n, float lat, float etpSeries[])