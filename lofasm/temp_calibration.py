#! /usr/bin/env python

from lofasm.bbx import bbx
from lofasm.time import UT_to_LST
from lofasm.galaxy_model import galaxyPower
from lofasm.station import LoFASM_Stations
from astropy.time import Time
from lofasm.parse_data import freq2bin
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime, timedelta
from glob import glob
from lofasm.filter import running_median
import os


def least_squares(xdata, ydata):    #function to find the least squares and returns the calibrated fit
    ''' using model Xi = G(Xi' + N)  I = 1 / G
    xdata being the lofasm data and ydata being the galaxy model'''
    xdata = np.array(xdata, dtype= np.float64)
    ydata = np.array(ydata, dtype= np.float64)

    m = len(xdata)
    Sxy = np.sum(xdata * ydata)
    Sxx = np.sum(xdata ** 2)
    Sx = xdata.sum()
    Sy = ydata.sum()


    #we concluded with this formula after solving the theres parameters to fit to...

    N = ((Sxx * Sy)-(Sxy * Sx))/ ((Sx ** 2) - (Sxx *  m))
    I = (Sy + (N * m)) / Sx
    # N = n
    #this Model

    calib = I * xdata - N

    return(calib, I, N)



if __name__ == '__main__':

    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('fnames',nargs='+')    #reading the given arguments into
    p.add_argument('freq',type = float)   #into variables
    p.add_argument('hor', type = int)
    args = p.parse_args()
    args.fnames.sort() #sort the files


    x = []

    start_hr_arr = np.zeros(len(args.fnames))

    fil = args.fnames[0]
    lf = bbx.LofasmFile(fil)
    st = int(lf.header['station'])


    # lsts, pows = galaxyPower.loadfile(st, args.freq, args.hor)

    for ii, f in enumerate(args.fnames):               #for loop to read through the data
        lf = bbx.LofasmFile(f)          #and append to an empty list
        lf.read_data()                  #the time in lst, the average power
        lofasm_pows = lf.data[:,freq2bin(args.freq)]    #and the given galaxy power
        avg = np.mean(lofasm_pows)                      #at that time

        #LOFASM 3 [
        # start = float(lf.header['start_mjd'])
        # start_postmjd = Time(start, format = 'mjd')
        # st = int(lf.header['station'])
        # lfstation = LoFASM_Stations[st]
        # lon = lfstation.lon * 180 / np.pi
        # start_sr = start_postmjd.sidereal_time('mean', lon)
        # start_hr = start_sr.hour
        # start_hr_arr[ii] = start_hr
        # #powmodel = galaxyPower.calculatepower(start_hr, lsts, pows)
        # ]

        # LOFASM 1 [
        start = lf.header['start_time']
        start_cv = Time(start, scale = 'utc')
        st = int(lf.header['station'])
        lfstation = LoFASM_Stations[st]
        lon = lfstation.lon * 180 / np.pi
        start_sr = start_cv.sidereal_time('mean', lon)
        start_hr = start_sr.hour
        lsts, pows = galaxyPower.loadfile(st, args.freq, args.hor)
        powmodel = galaxyPower.calculatepower(start_hr,lsts,pows)
        x.append((start_hr,avg,powmodel))
        # ]

        # x.append((start_hr,avg))

    # fout_name = 'xdata_{}'.format(args.freq)  #not to sure what we're using this for ...
    x = np.array(x)
    x = x.reshape(len(args.fnames),3)    #reshaping the appended list to be able to zip the data



    station = lf.header['station']   #this is for later naming the graphs
    start_date = start.split('T')[0]


    #fout_name = 'xdata_{}'.format(args.freq)  #used to read data to use again with running this whole script
    # x = np.array(x)
    # x = x.reshape(len(args.fnames),2)    #reshaping the appended list to be able to zip the data

    # gal = np.interp(start_hr_arr, lsts, pows)

    #this is for later naming the graphs
    # LOFASM 3 [
    # date_fmt = "%Y%m%d"
    # start_date = start_postmjd.datetime.strftime(date_fmt)
    # time_fmt = "%H%M%s"
    # start_time = start_postmjd.datetime.strftime(time_fmt)
    # ]

    t, avg, gal = zip(*x)

    calib, I, N = least_squares(avg,gal) #finding the parameters to calibrate the data
                                        #good to note that it is first being passed through a filters
                                        #and then running the calibration


    # calib_fixed_n, _, _ = least_squares_fixed_n(avg, gal)
    plt.figure()
    plt.title('Power vs LST ' + 'Station: ' + str(st) + ' Frequency: ' + str(args.freq) + 'Horizon: ' + str(args.hor))
    # plt.plot(t, 10*np.log10(avg), '.', label = 'LoFASM') #plots the raw lofasm data
    plt.plot(t, 10*np.log10(gal), '.', label = 'Model') #plots the galaxy model
    plt.plot(t, 10*np.log10(calib), '.', label = 'Calib I: {} N: {}'.format(I, N)) #plots the calibrated data



    avg_med = running_median(avg, 20) #hm seems like the first median filter in ln129 should not be there

    calib_filt, I, N = least_squares(avg_med,gal) #calibration with the median filter


    plt.plot(t, 10*np.log10(calib_filt), '.', label = 'Filtered Calib') #plot with filtered calibration


    plt.xlabel("LST (hrs)")
    plt.ylabel("Power (dB)")
    plt.legend()
    plt.savefig('Lofasm-Calib_' + str(st) + '_' + start_date + '_' + str(args.freq) + '_' + str(args.hor) + '.png')


    file = open("parameters.txt", 'a') #the next two lines are appending to file where the Station
                                       #frequency and time are recorded as well as the calibration parameters, I and N

    # file.write("Station: {}\tFrequency: {}\tHorizon: {}\nTime Stamp: {}\tI: {}\tN: {}\n".format(st, str(args.freq), str(args.hor), start_time, str(I), str(N)))

    file.close()
