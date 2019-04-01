import lofasm.calibrate.lofasmcal as lfc
import numpy as np
from datetime import datetime
from datetime import timedelta
from astropy.time import Time
from lofasm import parse_data as pdat
galaxyobj = lfc.galaxy()






lst_times = np.linspace(0,24,300)

lst_time = lst_times.tolist()
# galaxymodel = galaxyobj.galaxy_power_array(dates, freq, stationid, horizon_cutoff_alt)



# the output should be an array of data points with timestamp, and put that in a file
# every five minutes for 24 hours

hor = 0

st = 3


# while hor < 25:
    # 288: 5 minutes every hour for 24 hours
    # while st < 5: #loop through stations 1,2,3 and 4
    #     if z == 2: # skipping station2
    #         z = z + 1
    # df = 100./1024 #delta frequency
    # cur = df/2 #cursor for the average frequency in a bin
cur = 15.
df = 70. / 300.
# surya edit
df = 70./1024
# surya edit off
while cur < 86: # looping through range of frequencies 5-85 : UPDATE, changed to 101 to models at all frequencies across the 1024 spectra



    datafile = open("lofasm" + str(st) + "_" + str(cur) + "MHz_" + str(float(hor)) + "deg_new.dat", 'w')
    datafile.write("<Time Stamp>\t<Galaxy Power>\n" )


    galaxymodel = galaxyobj.galaxy_power_array(lst_time, cur, st, hor,lst=True)

    for x in range(len(lst_time)): #looping through times, to print it out

        datafile.write(str(lst_time[x]) + " " + str(galaxymodel[x]) + '\n')


    datafile.close()
    # print("freq: ", cur)
    cur = cur + df  #change between every frequency channel

        #z = z + 1
    # print("horizon: ", hor)
    # hor = hor + 5
