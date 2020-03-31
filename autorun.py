#!/usr/bin/python
# -*- coding:utf-8 -*-
import os
import sys

def run(cmd, file_name, times):
    times = int(times)
    pos_sum = 0
    best_pos = 15
    point_sum = 0
    for i in range(times):
        try:
            os.system(cmd)
        except Exception:
            print(f'round {i} passed, error occured')
            continue
        with open(file_name, 'r') as fin:
            result = fin.readlines()[-1].split(' ')
            pos = int(result[2])
            if pos < best_pos:
                best_pos = pos
            point = int(result[3])
        pos_sum += pos
        point_sum += point
        pos_average = pos_sum / (i+1)
        point_average = point_sum / (i+1)
        print(f'round {i}, ranked {pos}, point {point}')
        print(f'average rank {pos_average}, average point {point_average}, best rank {best_pos}')
    pos_average = pos_sum / times
    point_average = point_sum / times
    print(f'average rank {pos_average}, average point {point_average}, best rank {best_pos}')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f'usage: python autorun.py <cmd> <result file path> <times to simulate>')
        print(f'e.g. python autorun.py 2020_C ../xxx/result.txt, 1000')
    else:
        args = sys.argv
        run(args[1], args[2], args[3])
