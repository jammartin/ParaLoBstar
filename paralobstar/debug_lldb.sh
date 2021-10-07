#!/usr/bin/env bash

mpirun -np 2 xterm -e lldb bin/paralobstar -s initPipe.lldb
