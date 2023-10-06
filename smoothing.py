import cv2 as cv
img= cv.imread('cat1.jpg')
cv.imshow('Cat1',img)
average= cv.blur(img, (7,7)) #Averaging- here (7,7) is size of kernel averages all values around the center pixel
cv.imshow('Average Blur',average)
gauss= cv.GaussianBlur(img, (7,7),0)
cv.imshow('Gauss', gauss)
median= cv.medianBlur(img,3) #Median Blur
cv.imshow('Median Blur',median)
bilateral=cv.bilateralFilter(img,5,15,15)
cv.imshow('Bilateral',bilateral)
cv.waitKey(0)