import os
import ewok

__all__ = ["fc_file", "obs_file", "r2d2_obsfile", "r2d2_anfile"]


def fc_file(fcout, step):
    fc = {}
    fc['date'] = fcout['date']
    keys = [fcout['exp'], fcout['type'], fcout['date'], ewok.jediformat(step)]
    fname = '.'.join(keys)
    fc['filename'] = os.path.join(fcout['datadir'], fname)
    return fc


def obs_file(conf):
    obsfile = conf['obsdatain']
    return obsfile


def r2d2_obsfile(conf, date):
    sdate = ewok.jediformat(date)
    r2d2keys = ['l95', conf['source'], sdate, 'obt']
    r2d2file = '.'.join(r2d2keys)
    return r2d2file


def r2d2_anfile(conf, date):
    sdate = ewok.jediformat(date)
    r2d2keys = [conf['exp'], conf['type'], sdate, 'l95']
    r2d2file = '.'.join(r2d2keys)
    return r2d2file

