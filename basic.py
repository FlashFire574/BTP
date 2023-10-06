import cv2 as cv
import numpy as np
img= cv.imread('cat1.jpg')
 #Coverting to Grayscale
gray = cv.cvtColor(img,cv.COLOR_BGR2GRAY)
cv.imshow("Gray",img)

#Blur
blur=cv.GaussianBlur(img,(3,3), cv.BORDER_DEFAULT)
cv.imshow("Blur",img)

#grab edges of images
canny= cv.Canny(img,125,175)
cv.imshow("Canny",canny)

#resize
resized=cv.resize(img,(500,500))
cv.imshow("Resized",img)
cv.waitKey(0)
