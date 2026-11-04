"""
Scripts to load DAS HDF5 files created by Pascal Edme (e.g. after downsampling),
Derived from Sara and Daniel's original code to read Silixa iDASv2 format 
Pascal Edme, ETH Zürich
Last updated: Ju1ly 2023
"""

import os
import pathlib
import sys
import warnings

import h5py
import numpy as np
from datetime import datetime, timedelta
import glob
import h5py
import pathlib
import sys
import os

l_fields = []
l_attrs = []


def browse_file_attributes(name):
    l_fields.append(str(name))
    return None


def reset_attributes():
    l_fields.clear()
    l_attrs.clear()
    return None


def load_headers_only(file, verbose=False):
    """
    headers = load_das_h5.load_headers_only( file )

    :Return only the header, metadata, without the large data block
    :OUTPUTS:
    :headers -- dict containing many header information, such as:
    :   fs    -- sample rate in Hz
    :   dx    -- reported channel spacing in meters (may be off by -2%)
    :   lx    -- total fiber length
    :   nchan -- number of channels
    :   npts  -- number of samples
    :   t0    -- datetime object, start time of the file
    :   t1    -- datetime object, end time of the file
    :   d0    -- starting channel distance along fibre (usually the first ~200m are internal to the DAS)
    :   d1    -- end channel distance along fibre
    :   fm    -- fiber multiplier. Requested channel spacing of 1m is actually this distance
    :   unit  -- String of the units of measurment (i.e., strain-rate is (nm/m)/s * Hz/m)
    :   amp_scaling -- by default no scaling; i.e., 1.0
    :                  raw Silixa units are proportional to strain rate, but need scaling
    :                  to be accurate. This variable tracks what scaling value has been
    :                  applied. (see "convert=True" flag on the custom reader below)
    """
    # -- Reset containers for headers
    reset_attributes()

    with h5py.File(file, "r") as f:
        headers = dict()
        # headers["fs"] = int(f.attrs["Freq"])
        headers["fs"] = 1000 / float(f.attrs["step"])
        if "scaling" in f.attrs:
            headers["amp_scaling"] = f.attrs["scaling"]
        else:
            print("Amp scaling does not exist as attribute: set here to 1 !!!")
            headers["amp_scaling"] = 1

        if "gaugeLength" in f.attrs:
            headers["GaugeLength"] = f.attrs["gaugeLength"]
            print("GaugeLength attribute found (=", headers["GaugeLength"], "m)")

        ## can add some custom attributes here

        headers["fm"] = 1 #1.02095199
        headers["unit"] = "(m/m)/s"
        if verbose:
            print("step is", int(f.attrs["step"]))

        # # Check if "offset" is in the header_info array
        # headers["header_info"] = f.attrs["header_info"][0, :]
        # if "offset" in headers["header_info"]:
        #     # Find the indices where "offset" is present in the array
        #     offset_index = np.where(headers["header_info"] == "offset")[0]
        #     if len(offset_index) > 0:
        #         print("'offset' index in the header_info array:", offset_index)
        #     else:
        #         print("'offset' not found in the array.")
        # else:
        #     print("'offset' not found in the array.")

        # if "header_info" in f.attrs:
        #     headers["header_info"] = f.attrs["header_info"][0, :]
        #     print("Type of header_info", type(headers["header_info"]))
        #     print("shape of header_info", headers["header_info"].shape)
        #     print("header_info attribute found (=", headers["header_info"], "xxx)")
        #     if "offset" in headers["header_info"]:
        #         print("YYYY")
        #         position = headers["header_info"].index("offset")
        #         print("offset index", poisition)

        # -- Get dimensions from data array
        dset = f["Traces"]
        if verbose:
            print("data format is", dset.dtype)
        ddims = dset.shape
        headers["npts"] = ddims[0]
        headers["nchan"] = ddims[1]
        if verbose:
            print(
                "file size is",
                headers["npts"],
                " samples by",
                headers["nchan"],
                " channels with fs =",
                headers["fs"],
                "Hz i.e. duration=",
                headers["npts"] / headers["fs"],
                "sec",
            )

        # -- get headers dataset (at least the offset array here)
        # Check if "offset" is in the header_info array

        headers["header_info"] = f.attrs["header_info"][0, :]
        if "offset" in headers["header_info"]:
            # Find the indices where "offset" is present in the array
            offset_index = np.where(headers["header_info"] == "offset")[0]
            if verbose:
                if len(offset_index) > 0:
                    print("'offset' index in the header_info array:", offset_index)
                else:
                    print("'offset' not found in the array.")

            if "headers" in f:
                if verbose:
                    print("offset_index=", int(offset_index))
                dset = f["headers"][int(offset_index), :] / headers["amp_scaling"]
                # dset = f["headers"][0, :] / headers["amp_scaling"]
                if verbose:
                    print("Shape of header array dataset:", dset.shape)
                headers["d0"] = dset[0]  # first value
                headers["d1"] = dset[-1]  # last value
                headers["dx"] = dset[1] - dset[0]  # diff

            if verbose:
                print(
                    "Distance min=",
                    headers["d0"],
                    "Distance max=",
                    headers["d1"],
                    "dx=",
                    headers["dx"],
            )
        else:
            if verbose:
                print("'offset' not found in the header_info array attribute.")

        # # -- Get offset from Distance dataset
        # if "Distance" in f:
        #     dset = f["Distance"]
        #     headers["d0"] = dset[0]  # first value
        #     headers["d1"] = dset[-1]  # last value
        #     headers["dx"] = dset[1] - dset[0]  # diff
        #     print("Distance dataset exists.")
        # else:
        #     print("Distance dataset does not exists. try to read attribute offset")
        #     if "offset" in f.attrs:
        #         offset = f.attrs["offset"]
        #         headers["d0"] = offset[0]  # first value
        #         headers["d1"] = offset[-1]  # last value
        #         headers["dx"] = offset[1] - offset[0]  # diff
        #     else:
        #         print("Distance dataset does not exists nor attribute offset")
        #         print(">>> Cable distance set as Nchannel (dx=1)")
        #         headers["d0"] = 0
        #         headers["d1"] = headers["nchan"] - 1
        #         headers["dx"] = 1

        # headers["d0"] = headers["d0"] * headers["fm"]
        # headers["d1"] = headers["d1"] * headers["fm"]
        # headers["dx"] = headers["dx"] * headers["fm"]
        # print(
        #     "Distance min=",
        #     headers["d0"],
        #     "Distance max=",
        #     headers["d1"],
        #     "dx=",
        #     headers["dx"],
        # )



        # -- Infer times from the filename
        h5name = os.path.basename(file)
        headers["t0"] = datetime.strptime(h5name[-22:], "%Y%m%d_%H%M%S.%f.h5")
        headers["t1"] = headers["t0"] + timedelta(
            seconds=headers["npts"] / headers["fs"]
        )
        if verbose:
            print(
                "Loading headers from: {0}, start: {1}, end: {2}".format(
                    file, headers["t0"], headers["t1"]
                )
                )
        # if verbose:
        #     print(
        #         "Loading headers from: {0}, start: {1}, end: {2}".format(
        #             file, headers["t0"], headers["t1"]
        #         )
        #     )

    return headers


def load_file(file, convert=False, return_axis=True, verbose=False):
    """
    data, headers, axis = load_das_h5.load_file( file )

    :Load a single HDF5 file, based solely on the filename
    :OUTPUTS:
    :data    -- 2D numpy array [ num_samples, num_channels ]
    :headers -- dict of header information
    :axis    -- dict of constructed axis vectors:
    :            - tt         -- timesteps in seconds
    :            - date_times -- absolute datetime objects
    :            - dd         -- channel distances
    """

    headers = load_headers_only(file)
    # -- (Using this other function to load headers is cleaner, but currently
    # --   it means the file is opened and closed twice)

    with h5py.File(file, "r") as f:
        data = f["Traces"][:]
        data = data / headers["amp_scaling"]
        headers["unit"] = "(m/m)/s"
        if verbose:
            print(
                "Converted to strain rate (m/m)/s with amp_scaling = ",
                headers["amp_scaling"],
            )

    # # -- Convert everything
    # if convert == True:
    #     # -- DCB note: Silixa raw PRODML files report units as strain rate, even when
    #     # --  the following conversion has not yet been applied. Up to the user to mark
    #     # --  whether data has actually been converted yet or not...
    #     # -- We started writing a custom header defining any amplitude scaling that has been
    #     # --  applied (i.e., 1.0 if no scaling, some number otherwise)
    #     scale = True
    #     if "amp_scaling" in headers:
    #         # -- Check if amplitude has been scaled yet or is still 1.0
    #         # --  (close to 1.0, possible rounding errors with read/write)
    #         if np.abs(headers["amp_scaling"] - 1.0) > 0.0001:
    #             print(
    #                 'WARNING: flag "convert" is TRUE, but units are already scaled somehow'
    #             )
    #             print("   Doing nothing regarding conversion.")
    #             scale = False

    #     # -- DCB note: the sample rate (fs) used below is the ORIGINAL sample rate
    #     # --  at which data is recorded. If files have been downsampled, use the custom
    #     # --  header['fs_orig']
    #     if scale:
    #         # if "fs_orig" in headers.keys():
    #         #     fs = headers["fs_orig"]

    #         data = data / headers["amp_scaling"]
    #         # headers["amp_scaling"] = 116.0 / 8192.0 * fs / 10.0
    #         headers["unit"] = "(m/m)/s"
    #         print("Converted to strain rate!")
    #         # if verbose:
    #         #     print("Converted to strain rate!")

    ##########################
    # -- One could stop here...
    # -- But we'll compute vectors for channel spacing and timing
    if return_axis:
        axis = dict()

        # -- For the SWP Silixa iDAS, 1m channel spacing is actually ~1.02. See "FibreLengthMultiplier"
        # -- For the SWP Silixa iDAS, the first ~200m of fiber is internally inside the interrogator
        # --    See "StartDistance" for the exact value (which can change depending on settings)
        d0 = headers["d0"]
        d1 = headers["d1"]
        fm = headers["fm"]
        dx = headers["dx"]

        # dd = np.arange(d0, d1+dx*fm, dx*fm)
        # dd = np.arange(d0, d1, dx*fm)
        # -- NOTE: The end channel location calculated this way can differ from "d1" ("StopDistance") by
        # --  order 0.01m, especially when accumulated over a long (>30km) fibre. Rounding erorrs?
        # -- Possibly one would need to add/subtract one index to dd to get the dimensions correct.
        # -- Temporary solution? Add 1/2 a sample to the end target of np.arange, to make sure the final sample is reached
        # dd = np.arange(d0, d1 + dx * fm / 2, dx * fm)  # orig
        dd = np.arange(d0, d1 + dx * fm / 2, dx)  # pascal
        if len(dd) != np.shape(data)[1]:
            print(
                "WARNING. Data returned does not match channel spacing calculated from metadata (dd)"
            )
            print(
                "dd is {0} channels from {1} to {2}, spacing {3}".format(
                    len(dd), dd[0], dd[-1], dd[1] - dd[0]
                )
            )
            print("Reported metadata: d0: {0}, d1:{1}".format(d0, d1))
            print("data is {0} channels".format(np.shape(data)[1]))
        axis["dd"] = dd

        # -- Set up vectors of time samples
        fs = headers["fs"]
        t0 = headers["t0"]
        tt = np.arange(0, np.shape(data)[0] / fs, 1.0 / fs)
        date_times = [None] * len(tt)
        for i, t in enumerate(tt):
            # date_times.append((final_t0+t)._get_datetime())
            date_times[i] = t0 + timedelta(seconds=t)
        axis["tt"] = tt
        axis["date_times"] = date_times

    # print(file, t0, t1)
    if return_axis:
        return data, headers, axis
    else:
        return data, headers


def load_das_custom(
    t_start,
    t_end,
    d_start=0,
    d_end=0,
    convert=False,
    verbose=False,
    input_dir="./",
    return_axis=True,
    nth_channel=1,
):
    """
    data, heades, axis = load_das_custom(t_start, t_end, d_start=0, d_end=0, convert=False, verbose=False, input_dir='./')
    :Custom function to load files in a flexible way.
    :User defines a start time and end time, and data will be returned accordingly.
    :This is NOT sophisticated...
    : -data filenames must contain the date and time in a certain format
    : -filenames must be sortable according to date
    : -gaps in the data will not be handled correctly
    : -it *should* be able to handle requests that span multiple days, but no guarantees
    : -etc.
    :
    :It is intended to be faster than loading full blocks. By knowing the indices of data needed,
    : one can avoid loading the entire data block in HDF5.
    :
    :INPUTS:
    :t_start -- datetime object the request should start at
    :              (e.g. t_start = datetime.strptime('2021/10/12 09:24:30.0', '%Y/%m/%d %H:%M:%S.%f'))
    :t_end   -- datetime object end
    :d_start -- (optional) fibre distance at which to start returning data
    :d_end   -- (optional) fibre distance at which to stop
    :convert -- (optional) boolean to convert to strain rate if not already done
    :            WARNING: This requires knowing the sample rate of the raw data.
    :            If data had been downsampled but not converted, you will need to change "fs" in that conversion
    :verbose -- (optional) boolean to print more information about what is being loaded
    :input_dir -- (optional) string of directory in which to look for data
    :
    :OUTPUTS:
    :data    -- 2D numpy array [ num_samples, num_channels ]
    :headers -- dict of header information
    :axis    -- dict of constructed axis vectors:
    :            - tt         -- timesteps in seconds
    :            - date_times -- absolute datetime objects
    :            - dd         -- channel distances
    """

    ##############################################
    ## STEP 1: Determine possible files to load, based on directory and filename
    ## This assumes a standard filename contains, somewhere, the date and time
    ##  fomatted like "%Y%m%d_%H%M%S" (20211119_124730).
    ##############################################
    # -- All the files in this directory
    all_files = sorted(glob.glob(input_dir + "/*.h5"))
    file_gap = timedelta(minutes=5)  # -- If we have a gap of more than 5 minutes, we'll assume it's a new file

    # -- Check further, if we're requesting a span of more than one day
    # -- e.g., if t_end_dir is onto the next day or multiple days
    # -- We'll check a potential range of days. Specify the start day and end day in case they're different
    t_start_dir = (t_start - file_gap).strftime("/%Y_%m_%d/")
    t_start_day = datetime.strptime(t_start_dir, "/%Y_%m_%d/")
    t_end_dir = (t_end + file_gap).strftime("/%Y_%m_%d/")
    t_end_day = datetime.strptime(t_end_dir, "/%Y_%m_%d/")

    t_step_day = t_start_day
    while t_step_day <= t_end_day:
        t_step_dir = t_step_day.strftime("/%Y%m%d/") ### this is the subdirectory format
        # print(t_step_dir)
        if os.path.exists(input_dir + t_step_dir):
            all_files += sorted(glob.glob(input_dir + t_step_dir + "*.h5"))
        t_step_day += timedelta(days=1)

    all_files = sorted(all_files)

    # -- To avoid loading every header in the directory, do a quick search of candidates based on filename
    # -- We're only searching based broadly on the minutes listed, rather than trying to match filenames exactly
    consider_files = []
    t_step = t_start - file_gap

    while t_step <= t_end:
        ## glob seems slower than my for loop over the previously indexed/globbed array
        # consider_files += glob.glob("{0}/*{1}*.h5".format(input_dir,t_step.strftime('%Y%m%d_%H%M')))
        consider_string = "{0}".format(t_step.strftime("%Y%m%d_%H%M"))
        for file in all_files:
            if consider_string in file:
                consider_files += [file]
        t_step += timedelta(minutes=1)

    # -- If using longer (i.e., 10 min blocks), then we might miss a file this way
    # -- Kludgy fix? Repeat with hour resolution if needed
    if len(consider_files) == 0:
        if verbose:
            print(
                "Warning: no valid files found for consideration yet, trying to find anything within the day"
            )
        t_step = t_start
        while t_step < t_end:
            consider_string = "{0}".format(t_step.strftime("%Y%m%d"))
            print(consider_string)
            for file in all_files:
                if consider_string in file:
                    consider_files += [file]
            t_step += timedelta(hours=1)

    if len(consider_files) == 0:
        print("ERROR! No files found to try loading")
        print("  This could be an error in defining input_dir,")
        print("  or filenames don't match up with requested dates")
    ##############################################
    ## STEP 2: Skim headers from each file considered,
    ##  if timewindow is valid then load the data as we go
    ##############################################
    if verbose:
        print("----------------------------------------------")
        print(
            "----- Requested: {0} to {1}".format(
                t_start.strftime("%Y-%m-%d %H:%M:%S"),
                t_end.strftime("%Y-%m-%d %H:%M:%S"),
            )
        )
        print("----- Checking through files, loading as we go")
        print("----------------------------------------------")

    found_data_yet = False
    for filename in consider_files:
        # -- Reset containers for headers
        reset_attributes()

        # -- Open the headers of each file and look at the times
        headers = load_headers_only(filename, verbose=verbose)
        t0 = headers["t0"]
        t1 = headers["t1"]
        npts = headers["npts"]
        nchan = headers["nchan"]
        fs = headers["fs"]
        dx = headers["dx"]
        fm = headers["fm"]
        d0 = headers["d0"]
        d1 = headers["d1"]

        # -- Check if either t0 or t1 (or both) lies within the desired bounds of t_start and t_end
        # print(t_start.strftime('%Y-%m-%d %H:%M:%S.%f'))
        # print(t_end.strftime('%Y-%m-%d %H:%M:%S.%f'))
        # print(t0.strftime('%Y-%m-%d %H:%M:%S.%f'))
        # print(t1.strftime('%Y-%m-%d %H:%M:%S.%f'))
        ## Do we use any of the data in this file? Consider the 4 cases for requests
        ##   t_0                                               t_1
        ##    |***************file******************************|
        ##
        ## |---case1---|
        ## ts         te
        ##                                                    |----case2---|
        ##                                                    ts          te
        ##                      |---case3---|
        ##                      ts          te
        ##
        ## |-------------------------------case4----------------------------|
        ## ts                                                               te
        ##
        ##      # case 1 & 4                    # case 2 & 4             # case 3
        if (
            ((t_start < t0 and t0 < t_end) or (t_start < t1 and t1 < t_end))
            or (t0 < t_start and t_start < t1)
            or (t0 == t_start)
            or (t1 == t_end)
        ):
            if verbose:
                print("Use it!")

            # -- Define the time index from which to pull
            tt = np.arange(0, npts / fs, 1.0 / fs)
            # -- Initial values: full range
            i_pull_start = 0
            i_pull_end = npts - 1

            if t_start > t0:  # See if we should pull less on the front end
                t_rel = (t_start - t0).total_seconds()
                i_pull_start = np.argmin(np.abs(tt - t_rel))
                if verbose:
                    print("~~~ cut front ~~~~~~~~")
                    print(filename)
                    print("Requested start: {0}".format(t_start))
                    print("This file start: {0}".format(t0))
                    print("Starting at {0} seconds in".format(tt[i_pull_start]))

            if t1 > t_end:  # See if we should cut some off the end
                t_rel = (t_end - t0).total_seconds()
                i_pull_end = np.argmin(np.abs(tt - t_rel))
                if verbose:
                    print("~~~ cut end   ~~~~~~~~")
                    print(filename)
                    print("Requested end: {0}".format(t_end))
                    print("This file end: {0}".format(t1))
                    print("Cutting at {0} seconds in".format(tt[i_pull_end]))

            # -- Set up axis of distances / channels
            # dd = np.arange(d0, d1+dx*fm, dx*fm)
            # dd = np.arange(d0, d1, dx*fm)
            # -- NOTE: The end channel location calculated this way can differ from "d1" ("StopDistance") by
            # --  order 0.01m, especially when accumulated over a long (>30km) fibre. Rounding erorrs?
            # -- Possibly one would need to add/subtract one index to dd to get the dimensions correct.
            # -- Temporary solution? Add 1/2 a sample to the end target of np.arange, to make sure the final sample is reached
            dd = np.arange(d0, d1 + dx * fm / 2, dx * fm)

            # -- Open the file again for actual reading of data
            with h5py.File(filename, "r") as f:
                # -- Did the user specify any cutting along distance axis?
                if d_end > 0:
                    id1 = np.argmin(np.abs(dd - d_start))
                    id2 = np.argmin(np.abs(dd - d_end))
                    ###### CUSTOM
                    if nth_channel > 1:
                        d_pull = np.arange(id1, id2 + 1)
                        d_pull = d_pull[::nth_channel]
                        if verbose:
                            print("pulling every {0} traces".format(nth_channel))
                        data_tmp = f["Traces"][i_pull_start : i_pull_end + 1, d_pull]
                        dd = dd[d_pull]
                        headers["dx"] = headers["dx"] * nth_channel
                    ###### CUSTOM
                    else:
                        data_tmp = f["Traces"][
                            i_pull_start : i_pull_end + 1, id1 : id2 + 1
                        ]
                        dd = dd[id1 : id2 + 1]
                    if verbose:
                        print(
                            "Returning channels over distances: {0}  --  {1}".format(
                                dd[0], dd[-1]
                            )
                        )
                        print(
                            "data: {0} x {1},  dd: {2}".format(
                                np.shape(data_tmp)[0],
                                np.shape(data_tmp)[1],
                                np.shape(dd),
                            )
                        )
                # -- Otherwise just return all channels
                else:
                    if nth_channel > 1:
                        print(
                            "ERROR, nth_channel specified but d_start or d_end not specified."
                        )
                        print("not programmed to handle this situation.")
                        return
                    data_tmp = f["Traces"][i_pull_start : i_pull_end + 1, :]
                    if verbose:
                        print("Returning all channels")

            # See if we've defined data yet, otherwise concatenate with previous
            if found_data_yet == False:
                data = data_tmp
                final_t0 = t0 + timedelta(seconds=tt[i_pull_start])
                final_t1 = t0 + timedelta(seconds=tt[i_pull_end])
                found_data_yet = True
            else:
                data = np.concatenate((data, data_tmp))
                final_t1 = t0 + timedelta(seconds=tt[i_pull_end])

            # print(np.shape(data))
            # NOTE, growing large arrays like this is slow.
            # Better would be to determine all the useful times,
            #  init an emtpy array of the correct size,
            #  and then fill elements

    # print(final_t0)
    # print(final_t1)
    try:
        data
    except:
        print("ERROR! No data was loaded")

    ###################################
    ## STEP 3: Compute vectors for channel spacing and timing
    ## and update headers
    ##################################
    # -- Update headers to reflect the pulled data
    headers["d0"] = dd[0]
    headers["d1"] = dd[-1]
    headers["t0"] = final_t0
    headers["t1"] = final_t1
    headers["npts"] = np.shape(data)[0]
    headers["nchan"] = np.shape(data)[1]

    axis = dict()

    # -- dd was already created when looking up indices above
    # -- (No need to recalculate it at the end; this function assumes each file read has the same channels)
    if len(dd) != np.shape(data)[1]:
        print(
            "WARNING. Data returned does not match channel spacing calculated from metadata (dd)"
        )
        print(
            "dd is {0} channels from {1} to {2}, spacing {3}".format(
                len(dd), dd[0], dd[-1], dd[1] - dd[0]
            )
        )
        print("Reported metadata: d0: {0}, d1:{1}".format(d0, d1))
        print("data is {0} channels".format(np.shape(data)[1]))
    axis["dd"] = dd

    if return_axis:
        # -- Convert time-axis into numpy-happy date time objects
        # (also, not using the obspy UTCdatetime object)
        tt = np.arange(0, np.shape(data)[0] / fs, 1.0 / fs)
        date_times = [None] * len(tt)
        for i, t in enumerate(tt):
            # date_times.append((final_t0+t)._get_datetime())
            date_times[i] = final_t0 + timedelta(seconds=t)
        axis["tt"] = tt
        axis["date_times"] = date_times

    if verbose:
        print("----------------------------------------------")
        print(
            "Requested time chunk: {0}  --  {1}".format(
                t_start.strftime("%Y/%m/%d %H:%M:%S.%f"),
                t_end.strftime("%Y/%m/%d %H:%M:%S.%f"),
            )
        )
        print(
            "Returning time chunk: {0}  --  {1}".format(
                date_times[0].strftime("%Y/%m/%d %H:%M:%S.%f"),
                date_times[-1].strftime("%Y/%m/%d %H:%M:%S.%f"),
            )
        )
        print(
            "Data[ {0} samples, {1} channels ]".format(
                np.shape(data)[0], np.shape(data)[1]
            )
        )

    # -- Convert everything
    data = data / headers["amp_scaling"]
    if verbose:
        print(
            "in custom: Converted to strain rate (m/m/s) by division by amp_scaling =",
            headers["amp_scaling"],
        )

    # if convert == True:
    #     # -- DCB note: Silixa raw PRODML files report units as strain rate, even when
    #     # --  the following conversion has not yet been applied. Up to the user to mark
    #     # --  whether data has actually been converted yet or not...
    #     # -- We started writing a custom header defining any amplitude scaling that has been
    #     # --  applied (i.e., 1.0 if no scaling, some number otherwise)
    #     scale = True
    #     if "amp_scaling" in headers:
    #         # -- Check if amplitude has been scaled yet or is still 1.0
    #         # --  (close to 1.0, possible rounding errors with read/write)
    #         if np.abs(headers["amp_scaling"] - 1.0) > 0.0001:
    #             print(
    #                 'WARNING: flag "convert" is TRUE, but units are already scaled somehow'
    #             )
    #             print("   Doing nothing regarding conversion.")
    #             scale = False
    #     if scale == True:
    #         # -- DCB note: the sample rate (fs) used below is the ORIGINAL sample rate
    #         # --  at which data is recorded. If files have been downsampled, use the custom
    #         # --  header['fs_orig']
    #         if "fs_orig" in headers.keys():
    #             fs = headers["fs_orig"]

    #         data = 116.0 * data / 8192.0 * fs / 10.0
    #         headers["amp_scaling"] = 116.0 / 8192.0 * fs / 10.0
    #         headers["unit"] = "(nm/m)/s"
    #         if verbose:
    #             print("Converted to strain rate!")

    if return_axis:
        return data, headers, axis
    else:
        return data, headers
