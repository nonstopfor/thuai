#!/usr/bin/python
# -*- coding:utf-8 -*-
import matplotlib.pyplot as plt
import sys

if (len(sys.argv) != 3):
    print(f"usage: python {sys.argv[0]} <result file path> <player id>\ne.g.python {sys.argv[0]} round0_rank1.txt 15")
else:
    with open(sys.argv[1]) as fin:
        lines = fin.readlines()
    scores = []
    ranks = []
    for line in lines:
        words = line.split(' ')
        if (len(words) < 4):
            continue
        if words[0] == 'Player' and words[1] == sys.argv[2]:
            scores.append(int(words[3]))
            ranks.append(int(words[2]))
    x = [i+1 for i in range(len(scores))]
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.plot(scores)
    ax2 = ax1.twinx();
    plt.plot(ranks, 'r')
    plt.savefig(f'ranks_{sys.argv[2]}')
