# (C) Copyright 2020-2021 UCAR
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import os
import yamltools
import ewok.tasks.forecast as generic

# Specific task inherits from generic but this is not mandatory
class forecastQG(generic.forecast):

    def setup(self, config, execs, fix, ic):
        print("Showing off with a specific forecast for QG")

        # Setting generic defaults (but task could set its own)
        generic.forecast.setup(self, config, execs, fix, ic)

        # Using specific script (generic one could be used as well)
        self.command = os.path.join(config['model_path'],"tasks/qg-run.sh")

