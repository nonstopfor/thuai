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
    debuglog = 'debuglog'
    mkdir(debuglog)
    players_pos = [0 for i in range(16)]
    players_point = list(players_pos)
    players_pos_sum = list(players_pos)
    players_point_sum = list(players_pos)
    players_best_pos = [20 for i in range(16)]
    for i in range(times):
        try:
            os.system(cmd+" >> "+ debuglog+f"/round{i}_debug.txt")
        except Exception as e:
            output(f'round {i} passed, error occured: {str(e)}', record)
            continue
        with open(file_name, 'r') as fin:
            lines = fin.readlines()
            for j in range(len(players_pos)):
                players_pos[j] = int(lines[j-16].split(' ')[2])
                players_point[j] = int(lines[j-16].split(' ')[3])
                players_pos_sum[j] += players_pos[j]
                players_point_sum[j] += players_point[j]
                if players_pos[j] < players_best_pos[j]:
                    players_best_pos[j] = players_pos[j]
                pos_average = players_pos_sum[j] / (i+1)
                point_average = players_point_sum[j] / (i+1)
                best_pos = players_best_pos[j]
                if (j == 15 and ((pos >= pos_average and pos >= 6) or pos == 1))
                    with open(directory+f"/round{i}_rank{pos}.txt", 'w') as fout:
                        fout.writelines(lines)
                output(f'player {j}: round {i}, ranked {players_pos[j]}, point {players_point[j]}', record)
                output(f'\taverage rank {pos_average}, average point {point_average}, best rank {best_pos}',\
                record)
    record.close()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f'usage: python autorun.py <cmd> <result file path> <times to simulate>')
        print(f'e.g. python autorun.py 2020_C ../xxx/result.txt, 1000')
    else:
        args = sys.argv
        run(args[1], args[2], args[3])
