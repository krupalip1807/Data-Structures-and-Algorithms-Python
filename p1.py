#!/usr/bin/python3

import sys, random, argparse, time, datetime
import p1_algs


verbose = '''verbose levels: 0 - numeric answer only, 1 - include alignment, 2 - include DP table / lists used in computation'''

# set up the command-line arguments
parser = argparse.ArgumentParser(description = 'p1 - string algorithms toolkit')
group_x = parser.add_mutually_exclusive_group(required=True)
group_x.add_argument('--x', default = None, type=str, help='first string')
group_x.add_argument('--x_file', type=argparse.FileType('r'), help='file to read to use as first string')
group_y = parser.add_mutually_exclusive_group(required=True)
group_y.add_argument('--y_file', type=argparse.FileType('r'), help='second string')
group_y.add_argument('--y', default = None, type=str, help='file to read as second string')
parser.add_argument('--alg', default = '', choices= p1_algs.algorithms.keys(),
                    type=str, help=p1_algs.help, required=True)
parser.add_argument('--verbose', default=0, choices=[0,1,2],
                    type=int, help=verbose, required=False)
parser.add_argument('--all_matches', choices=['yes', 'no'], default='yes',
                    type=str, help='output all matches (yes) or only first (no) for string matching algorithms')
parser.add_argument('--reverse_y', choices=['yes', 'no'], default='no',
                    type=str, help='reverse y string before running algorithm')
parser.add_argument('--strip_space', default='yes', choices=['yes', 'no'],
                    type=str, help='remove space from strings before matching/alignment')
parser.add_argument('--comment_char', default='>', type=str, help='ignore lines starting with this char')

# parse the command-line arguments
args = parser.parse_args()

# strings will be X and Y
if args.x: X = args.x
else: X = args.x_file.read()

if args.y: Y = args.y
else: Y = args.y_file.read()

# modifications of X and Y

# comment char
if args.comment_char != '':
    def remove_comment(X, ch):
        X_list = X.split('\n')
        X = ''
        for i in range(0, len(X_list)):
            line = X_list[i]
            if len(line) > 0 and line[0] == ch: continue
            X = X + line
            if i < len(X_list)-1: X = X + '\n'
        return X
    X = remove_comment(X, args.comment_char)
    Y = remove_comment(Y, args.comment_char)

# strip space
if args.strip_space == 'yes':
    def rid_space(X):
        X = X.replace(' ','').replace('\t','').replace('\n','').replace('\r','')
        return X
    X = rid_space(X)
    Y = rid_space(Y)
    
# reverse
if args.reverse_y == 'yes':
    Y = Y[::-1]

# and call the appropriate algorithm based on the choice in algorithm
p1_algs.algorithms[args.alg](X, Y, args)
