#
# ====================================   SET_OTIME_FROM_FILENAME()  =======================================
#
def set_otime_from_filename(self, offset=0., prefix=None):
    """
    ##Description
    Set data header origin time from the filename assuming Febus convention on the filename
    ex: SR_2021-08-26_14-32-39_UTC.h5.
    This call affect the origin time to file header and is propagated to all subsequent A1File.read() call
    This does not affect the values of the <time> vector which always refer to the origin (ie starting) time
    ##Input
        offset = a time offset to add/substract from the filename information
        prefix = a prefix ending by "_" in case the filename do not follow Febus convention (SR_, RAW_, S_, ...)
    """
    from datetime import timezone, datetime
    from ._a1das_exception import WrongValueError

    name = self.file_header.filename
    # expected to be of the form SR_2021-08-26_14-32-39_UTC.h5
    if prefix is None:
        start = name.find("_")
        ofs = 1
    else:
        start = name.find(prefix)
        ofs = len(prefix)

    if start == -1:
        raise WrongValueError('could not set origin time from filename, check filename and/or prefix')

    end = name.find("_UTC")
    s = name[start + ofs:end]
    # convert from string
    try:
        d = datetime.strptime(s, "%Y-%m-%d_%H-%M-%S")
    except:
        raise ValueError('wrong date-time format, check file prefix')

    # convert to UTC
    dd = datetime(d.year, d.month, d.day, d.hour, d.minute, d.second, tzinfo=timezone.utc)
    # convert to Posix
    self.data_header.set_item(otime=dd.timestamp() + offset)


    #
    # ====================================   SET_OTIME_FROM_FILENAME()  =======================================
    #
    def set_otime_from_filename(self, offset=0., prefix=None):
        """
        ##Description
        Set data header origin time from the filename assuming Febus convention on the filename
        ex: SR_2021-08-26_14-32-39_UTC.h5.
        This call affect the origin time to file header and is propagated to all subsequent A1File.read() call
        This does not affect the values of the <time> vector which always refer to the origin (ie starting) time

        ##Input
            offset = a time offset to add/substract from the filename information
            prefix = a prefix ending by "_" in case the filename do not follow Febus convention (SR_, RAW_, S_, ...)
        """
        from datetime import timezone, datetime
        from ._a1das_exception import WrongValueError

        if 'filename' in self.data_header._header.keys():
            name = self.data_header['filename']
            # expected to be of the form SR_2021-08-26_14-32-39_UTC.h5
            if prefix is None:
                start = name.find("_")
                ofs = 1
            else:
                start = name.find(prefix)
                ofs = len(prefix)

            if start == -1:
                raise WrongValueError('could not set origin time from filename, check filename and/or prefix')

            end = name.find("_UTC")
            s = name[start+ofs:end]
            # convert from string
            try:
                d = datetime.strptime(s, "%Y-%m-%d_%H-%M-%S")
            except:
                raise ValueError('wrong date-time format, check file prefix')

            # convert to UTC
            dd = datetime(d.year,d.month,d.day,d.hour,d.minute,d.second, tzinfo=timezone.utc)
            # convert to Posix
            self.data_header.set_item(otime=dd.timestamp()+offset)
        else:
            raise WrongValueError('<filename> is not defined in data header')
