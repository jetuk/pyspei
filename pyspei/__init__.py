from ._spei import thornthwaite
import numpy as np

def read_txt(filename, delim=';', dtype=np.float64):
    """

    tampa
    27.96
    1900;01
    12
    110.70;14.30
    105.40;15.40
    186.90;18.00
    95.00;21.41
    """

    with open(filename, 'r') as fh:

        # Headers
        name = fh.readline().strip()
        latitude = float(fh.readline().strip())
        year, mth = map(int, fh.readline().strip().split(delim))
        seasonality = int(fh.readline().strip())

        precip, temp = [], []
        # remaining content is data
        for row in fh.readlines():
            p, t = map(float, row.strip().split(delim))

            precip.append(p)
            temp.append(t)

    return name, latitude, (year, mth), seasonality, np.array(precip, dtype=dtype), np.array(temp, dtype=dtype)


def spei(precip, interval=12, temp=None, latitude=None, seasonality=12, fit_only=False, return_fit_parameters=False,
         beta=None, logLogisticParams=None):
    from ._spei import spei as c_spei, fit_loglogistic
    pet = np.zeros_like(precip)
    # Compute PET using Thornthwaite equation
    if temp is not None and latitude is not None:
        thornthwaite(temp, latitude, pet)

    her = (precip - pet).astype(np.float64)

    if fit_only:
        # Only get and retun the fit parameters
        beta, logLogisticParams = fit_loglogistic(her, interval, seasonality)
        return np.array(beta), np.array(logLogisticParams)
    else:
        # Calculate SPEI
        spei_data = np.empty(precip.shape[0] - interval + 1, dtype=np.float64)
        beta, logLogisticParams = c_spei(her, interval, seasonality, spei_data,
                                         beta=beta, logLogisticParams=logLogisticParams)

    if return_fit_parameters:
        return spei_data, np.array(beta), np.array(logLogisticParams)
    else:
        return spei_data


def write_spei_txt(filename, spei, name, latitude, year, month, interval, delim=';'):

    with open(filename, 'w') as fh:

        fh.write('\n'.join([
            name, '{:.2f}'.format(latitude), '{:04d}{}{:02d}'.format(year, delim, month), '{:d}\n'.format(interval)
        ]))

        fh.writelines(
            ['{:.6f}\n'.format(v) for v in spei]
        )

