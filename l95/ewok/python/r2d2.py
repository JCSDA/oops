import os
import ewok

__all__ = ["obsfile", "fbfile", "anfile"]

def obsfile(conf, date):
    sdate = ewok.jediformat(date)
    r2d2keys = [conf['_source'], 'l95', sdate, 'obt']
    r2d2file = '.'.join(r2d2keys)
    return r2d2file

def fbfile(conf, date):
    sdate = ewok.jediformat(date)
    r2d2keys = [conf['experiment']['expid'], 'l95', 'fb', sdate, 'obt']
    r2d2file = '.'.join(r2d2keys)
    return r2d2file

def anfile(conf, date):
    sdate = ewok.jediformat(date)
    r2d2keys = [conf['expid'], 'l95', 'an', sdate, 'l95']
    r2d2file = '.'.join(r2d2keys)
    return r2d2file

def fcfile(conf, date, step):
    sdate = ewok.jediformat(date)
    sstep = ewok.jediformat(step)
    r2d2keys = [conf['expid'], 'l95', 'fc', sdate, sstep, 'l95']
    r2d2file = '.'.join(r2d2keys)
    return r2d2file

