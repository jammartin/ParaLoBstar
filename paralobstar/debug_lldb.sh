#!/usr/bin/env bash

mpirun -np 2 xterm -e lldb bin/runner -s initPipe.lldb
