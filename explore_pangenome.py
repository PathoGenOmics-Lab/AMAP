#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Explore pangenome')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input JSON file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output JSON file')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    

if __name__ == '__main__':
    main()