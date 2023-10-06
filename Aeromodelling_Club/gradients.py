import cv2 as cv
import numpy as np
img= cv.imread('cat1.jpg')
cv.imshow('Cat1',img)

gray=cv.cvtColor(img,cv.COLOR_BGR2GRAY)
cv.imshow('Gray',gray)

lap=cv.Laplacian(gray,cv.CV_64F) #image gradients a re transition of intensity black to white is + and white to dark is -
lap=np.uint8(np.absolute(lap))
cv.imshow('Lapaclain',lap)


sobelx=cv.Sobel(gray,cv.CV_64F,1,0)
sobely=cv.Sobel(gray,cv.CV_64F,0,1)
cv.imshow('Sobel X',sobelx)
cv.imshow('Sobel Y',sobely)
combined_sobel=cv.bitwise_or(sobelx,sobely)
cv.imshow('Combined Sobel',combined_sobel)

canny= cv.Canny(gray,150,175)
cv.imshow('Canny',canny)

cv.waitKey(0)