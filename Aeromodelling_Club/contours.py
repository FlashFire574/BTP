import cv2 as cv
import numpy as np
img= cv.imread('cat1.jpg')
cv.imshow('Cat1',img)

gray= cv.cvtColor(img,cv.COLOR_BGR2GRAY)
cv.imshow("Gray",gray)

blank = np.zeros(img.shape[:2], dtype="uint8")
cv.imshow("Blank", blank)



#blur = cv.GaussianBlur(gray,(5,5), cv.BORDER_DEFAULT)
#cv.imshow("Blur", blur)
#canny=cv.Canny(blur,125,175)
#cv.imshow("Canny",canny)

ret, thresh = cv.threshold(gray, 125, 255, cv.THRESH_BINARY) # binarizes the image
cv.imshow("Thresh", thresh)

contours, heirarchies = cv.findContours(thresh, cv.RETR_LIST , cv.CHAIN_APPROX_SIMPLE)
cv.drawContours(blank, contours, -1, (0,0,255), 2)
cv.imshow("Contours Drawn", blank)
cv.waitKey(0)
