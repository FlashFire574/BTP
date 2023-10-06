import cv2 as cv
import numpy as np

blank = np.zeros((400,400), dtype='uint8')
rectangle= cv.rectangle(blank.copy(),(30,30),(370,370),255, -1 ) # -1 is for filling
circle= cv.circle(blank.copy(), (200,200), 200, 255, -1) #-1 is for filling

cv.imshow('Rectangle', rectangle)
cv.imshow('Circle', circle)

bitwiseAND=cv.bitwise_and(rectangle, circle)
cv.imshow('Bitwise AND',bitwiseAND)

bitwiseOR=cv.bitwise_or(rectangle,circle)
cv.imshow('Bitwise OR', bitwiseOR)

bitwiseXOR=cv.bitwise_xor(rectangle, circle)
cv.imshow('Bitwise XOR',bitwiseXOR)
cv.waitKey(0)