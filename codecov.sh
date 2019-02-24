#!/bin/bash

CODECOV_TOKEN=$CODECOV_MODERNROBOTICSBOOK julia -e 'import ModernRoboticsBook; cd(joinpath(dirname(pathof(ModernRoboticsBook)), "..")); using Coverage; Codecov.submit_local(Codecov.process_folder())'
