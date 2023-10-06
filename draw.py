import cv2 as cv
import numpy as np
img= cv.imread('cat1.jpg')
cv.imshow('Cat1',img)
blank= np.zeros((500,500,3), dtype='uint8')
cv.imshow('Blank',blank)
#blank[:]= 0,255,0 #image colured green order is BGR
#cv.imshow('Green',blank)
#blank[200:300, 300:400]= 0,0,255
#cv.imshow('Red Square',blank)
cv.rectangle(blank,(0,0),(250,250),(0,255,0),thickness=10)
cv.waitKey(0)

