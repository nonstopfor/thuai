#!/usr/bin/python
# -*- coding:utf-8 -*-
import os
import sys

def run(cmd, file_name, times):
    times = int(times)
    pos_sum = 0
    point_sum = 0
    for i in range(times):
        os.system(cmd)
        with open(file_name, 'r') as fin:
            result = fin.readlines()[-1].split(' ')
            pos = int(result[2])
            point = int(result[3])
        pos_sum += pos
        point_sum += point
        print(f'round {i}, ranked {pos}, point {point}')
    pos_average = pos_sum / times
    point_average = point_sum / times
    print(f'average rank {pos_average}, average point {point_average}')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f'usage: python autorun.py <cmd> <result file path> <times to simulate>')
        print(f'e.g. python autorun.py 2020_C ../xxx/result.txt, 1000')
    else:
        args = sys.argv
        run(args[1], args[2], args[3])
