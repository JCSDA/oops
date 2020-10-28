import os
import yaml

class getBackground(Task):

  #  Need to register in task factory

  def setup(self, config, bg, obslist):
    bgconf = {}
    bgconf['date'] = '{{current_cycle}}'
    bgfile = os.path.join(config['datadir'],
                          '{{current_cycle}}',
                          'forecast.{{previous_cycle}}.{{window_offset}}')
    bgconf['filename'] = bgfile

    fnameout = config['workdir']/ewok/background.yaml
    yaml.dump(bgconf, fnameout)

    self.command = 'path/to/actual/runtime/qg/getbg.py'
    self.yamlarg = fnameout

    return bgconf
