import os
import ewok

__all__ = ["obsfile", "fbfile", "anfile"]

def obsfile(conf, date):
    sdate = ewok.jediformat(date)
    r2d2keys = [conf['_source'], 'qg', sdate, conf['obs type'], 'obs', 'nc']
    r2d2file = '.'.join(r2d2keys)
    return r2d2file

def fbfile(conf, date):
    sdate = ewok.jediformat(date)
    r2d2keys = [conf['experiment']['expid'], 'qg', 'fb', conf['obs']['_obsname'], sdate, 'obt']
    r2d2file = '.'.join(r2d2keys)
    return r2d2file

def anfile(conf, date):
    sdate = ewok.jediformat(date)
    r2d2keys = [conf['expid'], 'qg', 'an', sdate, 'nc']
    r2d2file = '.'.join(r2d2keys)
    return r2d2file

def fcfile(conf, date, step):
    sdate = ewok.jediformat(date)
    sstep = ewok.jediformat(step)
    r2d2keys = [conf['expid'], 'qg', 'fc', sdate, sstep, 'nc']
    r2d2file = '.'.join(r2d2keys)
    return r2d2file

