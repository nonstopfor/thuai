#!/usr/bin/python
# -*- coding:utf-8 -*-
import os
import sys

def mkdir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)

def output(s, fout):
    print(s)
    fout.write(s+'\n')
    fout.flush()

def run(cmd, file_name, times):
    times = int(times)
    pos_sum = 0
    best_pos = 15
    point_sum = 0
    record = open("./log.txt", 'w')
    directory = 'records'
    mkdir(directory)
    for i in range(times):
        try:
            os.system(cmd)
        except Exception:
            output(f'round {i} passed, error occured', record)
            continue
        with open(file_name, 'r') as fin:
            lines = fin.readlines()
            result = lines[-1].split(' ')
            pos = int(result[2])
            if pos < best_pos:
                best_pos = pos
            point = int(result[3])
        pos_sum += pos
        point_sum += point
        pos_average = pos_sum / (i+1)
        point_average = point_sum / (i+1)
        with open(directory+f"/round{i}_rank{pos}.txt", 'w') as fout:
            fout.writelines(lines)
        output(f'round {i}, ranked {pos}, point {point}', record)
        output(f'average rank {pos_average}, average point {point_average}, best rank {best_pos}',\
                record)
    pos_average = pos_sum / times
    point_average = point_sum / times
    output(f'final average rank {pos_average}, average point {point_average}, best rank {best_pos}',\
            record)
    record.close()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f'usage: python autorun.py <cmd> <result file path> <times to simulate>')
        print(f'e.g. python autorun.py 2020_C ../xxx/result.txt, 1000')
    else:
        args = sys.argv
        run(args[1], args[2], args[3])
