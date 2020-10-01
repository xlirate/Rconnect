#!/usr/bin/python3

#import struct

import secrets

width = 16
height = 16

block = list([list([secrets.randbelow(256) for _ in range(width)]) for _ in range(height)])

#print(block)

def int_to_bytes(n, b):
    arr = bytearray(b)
    for i in range(b):
        #arr[i] = (n//((2**8)**(b-(i+1))))%(2**8)
        arr[i] = (n//((2**8)**(i)))%(2**8)
    return arr





with open("./block.bmp", "wb") as bl:
    #Windows Structure: BITMAPFILEHEADER
    bl.write(b'BM')                   # Magic number
    bl.write(int_to_bytes(14+40+height*(width+width%4),  4))     # size, replace later
    bl.write(int_to_bytes(0,      2)) # unused
    bl.write(int_to_bytes(0,      2)) # unused
    bl.write(int_to_bytes(0x38,   4)) # offset to the start of the bitmap data

    #Windows Structure: BITMAPINFOHEADER
    bl.write(int_to_bytes(40,     4)) #header size
    bl.write(int_to_bytes(width,  4)) #width
    bl.write(int_to_bytes(height, 4)) #height
    bl.write(int_to_bytes(1,      2)) #magic number
    bl.write(int_to_bytes(8,      2)) #bits per pixel
    bl.write(int_to_bytes(0,      4)) #compression
    bl.write(int_to_bytes(0,      4)) #compressed size, 0 if not compressed
    bl.write(int_to_bytes(0,      4)) #pixels/meter horisontally
    bl.write(int_to_bytes(0,      4)) #pixels/meter vertically
    bl.write(int_to_bytes(0,      4)) #colours used,
    bl.write(int_to_bytes(0,      4)) #number of important colours, or 0 for all colours

    bl.write(bytes(2)) #padding

    for row in block:
        bl.write(bytes(row))
        bl.write(bytes(width%4))


