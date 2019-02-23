#!/bin/bash

CODECOV_TOKEN=$CODECOV_MODERNROBOTICS julia -e 'import ModernRobotics; cd(joinpath(dirname(pathof(ModernRobotics)), "..")); using Coverage; Codecov.submit_local(Codecov.process_folder())'
