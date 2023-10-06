import cv2 as cv
img= cv.imread('cat1.jpg')
cv.imshow('Cat1',img)
gray= cv.cvtColor(img,cv.COLOR_BGR2GRAY)
cv.imshow("Gray",gray)

hsv = cv.cvtColor(img,cv.COLOR_BGR2HSV)
cv.imshow("HSV",hsv)
cv.waitKey(0)