import os
from solo.date import JediDate, DateIncrement

__all__ = ["fc_file", "obs_file", "r2d2_obsfile", "r2d2_anfile"]


def fc_file(fcout, step):
    fc = {}
    fc['date'] = fcout['date']
    step = DateIncrement(duration=step)
    keys = [fcout['exp'], fcout['type'], fcout['date'], str(step)]
    fname = '.'.join(keys)
    fc['filename'] = os.path.join(fcout['datadir'], fname)
    return fc


def obs_file(conf):
    obsfile = conf['obsdatain']
    return obsfile


def r2d2_obsfile(conf, date):
    sdate = JediDate(date)
    r2d2keys = ['l95', conf['source'], str(sdate), 'obt']
    r2d2file = '.'.join(r2d2keys)
    return r2d2file


def r2d2_anfile(conf, date):
    sdate = JediDate(date)
    r2d2keys = [conf['exp'], conf['type'], str(sdate), 'l95']
    r2d2file = '.'.join(r2d2keys)
    return r2d2file

