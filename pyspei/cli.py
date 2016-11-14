# Command line interface for pyspei

from pyspei import read_txt, write_spei_txt, spei

def main(interval, input_filename, output_filename):

    name, lat, (year, mth), seasonality, precip, temp = read_txt(input_filename)

    print(precip.shape, interval)
    spei_data = spei(precip, interval=interval, temp=temp, latitude=lat)

    first_spei_year = year + (interval - 1) // 12
    first_spei_mth = mth + interval - 1
    while first_spei_mth > 12:
        first_spei_mth -= 12

    write_spei_txt(output_filename, spei_data, name, lat, first_spei_year, first_spei_mth, interval)


def cli():
    import argparse

    parser = argparse.ArgumentParser(description='Calculate SPEI from precipitation and temperature.')

    parser.add_argument('interval', type=int, help='Interval in months over which SPEI is calculated.')
    parser.add_argument('input_file', type=str, help='Input txt file containing precipitation and temperature data.')
    parser.add_argument('output_file', type=str, help='Output txt file to save SPEI.')

    args = parser.parse_args()

    main(args.interval, args.input_file, args.output_file)
