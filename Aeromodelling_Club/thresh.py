import cv2 as cv
img= cv.imread('cat1.jpg')
cv.imshow('Cat1',img)

gray=cv.cvtColor(img,cv.COLOR_BGR2GRAY)
cv.imshow('Gray',gray)

 #Simple Threshold
threshold,thresh=cv.threshold(gray,150,255,cv.THRESH_BINARY_INV)
cv.imshow('Simple Thresholded', thresh)

threshold,thresh=cv.threshold(gray,150,255,cv.THRESH_BINARY_INV)
cv.imshow('Simple Thresholded Inverse', thresh)

#Adaptive Threshold
adaptive_thresh=cv.adaptiveThreshold(gray,255,cv.ADAPTIVE_THRESH_MEAN_C,cv.THRESH_BINARY, 11,3) #blocksize is size of kernel required to take average on
cv.imshow('Adaptive Thresh',adaptive_thresh)
cv.waitKey(0)